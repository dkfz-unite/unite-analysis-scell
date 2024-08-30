# Unite Single Cell RNA Analysis Service

## General
Single cell RNA analysis service wrapperd with web API.


## Configuration
To configure the application, change environment variables as required in [commands](https://github.com/dkfz-unite/unite-commands/blob/main/README.md#configuration) web service:
- `UNITE_COMMAND` - command to run the analysis package (`python`).
- `UNITE_COMMAND_ARGUMENS` - command arguments (`app.py {data}/{proc}`).
- `UNITE_SOURCE_PATH` - location of the source code in docker container (`/src`).
- `UNITE_DATA_PATH` - location of the data in docker container (`/mnt/data`).
- `UNITE_LIMIT` - maximum number of concurrent jobs (`1` - process is heavy and uses a lot of CPU).


## Installation

### Docker Compose
The easiest way to install the application is to use docker-compose:
- Environment configuration and installation scripts: https://github.com/dkfz-unite/unite-environment
- Single cell RNA analysis service configuration and installation scripts: https://github.com/dkfz-unite/unite-environment/tree/main/applications/unite-analysis-sc

### Docker
[Dockerfile](Dockerfile) is used to build an image of the application.
To build an image run the following command:
```
docker build -t unite.analysis.sc:latest .
```

All application components should run in the same docker network.
To create common docker network if not yet available run the following command:
```bash
docker network create unite
```

To run application in docker run the following command:
```bash
docker run \
--name unite.analysis.sc \
--restart unless-stopped \
--net unite \
--net-alias sc.analysis.unite.net \
-p 127.0.0.1:5302:80 \
-e ASPNETCORE_ENVIRONMENT=Release \
-v ./data:/mnt/data:rw \
-d \
unite.analysis.sc:latest
```


## Usage

### Prepare The Data
Place samples data, metadata and processing options files to `{proc}` subdirectory of the `./data` directory on the host machine.
```txt
./data
└── {proc}
    ├── sample-(n)
        ├── features.tsv.gz
        ├── barcodes.tsv.gz
        └── matrix.mtx.gz
    ├── matadata.tsv
    └── options.json 
```
Where:
- `sample-(n)` - directories with single cell RNA data for each sample.
- `metadata.tsv` - metadata file where key column is `sample_id` and contains metadata for all required samples:
  ```tsv
  sample_id diagnosis sex age
  sample-1  glioblastoma  Male  45
  sample-2  glioblastoma  Female  55
  sample-n  glioblastoma  Male  65
  ```   
- `options.json` - processing options file with the following structure:
  ```jsonc
  {
    "qc": true, // Calculate quality control metrix.
    "sparce": true, // Make data sparce (highly recommended to set to 'true').
    "pp": "default", // Preprocessing option: default|seurat|zheng17.
    "pca": true, // Perform principal component analysis (required for neighbours calculation).
    "neighbours": true, // Perform neighbours calculation (required for clustering).
    "clustering": "louvain", // Clustering method: louvain|leiden.
    "embedding": "umap" // Embedding method: umap|tsne. 
  }
  ```

### Run The Analysis
Send a POST request to the `localhost:5302/api/run?key={key}` endpoint, where `key` is the process key and the name of the corresponding process directory.

This will invoke the command `python` with the arguments `app.py {data}/{proc}` where:
- All entries of `{proc}` will be replaced with the process `key`, which is the name of the corresponding process directory.
- All entries of `{data}` will be replaced with the path to the data location in docker container (In the example `./data` on the host machine will be mounted to `/mnt/data` in container).

### Analysis
Analysis will perform the following steps:
- Read the data for all samples from corresponding directories.
- Read the metadata from `metadata.tsv` file and enrich samples data with related metadata.
- Read the processing options from `options.json` file.
- Perform the analysis according to the processing options.
- Write resulting dataset to `data.h5ad` file.
