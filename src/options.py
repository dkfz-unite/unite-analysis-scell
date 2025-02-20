# class or struct with list of processing options
# qc: bool
# sparce: bool
# pp: str (default, 'seurat', 'zheng17')
# genes: int - required for default pp - minimum number of genes expressed in a cell
# cells: int - required for default pp - minimum number of cells with expressed gene
# annotate: bool - pefrom annotation
# model: str - required for annotation - model name
# pca: bool - required for neighbors
# neighbors: bool - required for clustering
# clustering: ('louvain', 'leiden')
# embedding: (['umap', 'tsne'])

import json

class Options:
    def __init__(
            self, 
            qc: bool = False,
            sparse: bool = True,
            # pp: str = 'default' | 'seurat' | 'zheng17',
            pp: str = 'default',
            genes: int = 5,
            cells: int = 25,
            annotate: bool = False,
            model: str = None,
            pca: bool = True,
            neighbors: bool = True,
            # clustering: str = 'louvain'|'leiden',
            clustering: str = 'louvain',
            # embedding: list = ['umap', 'tsne']):
            embedding: list = ['umap']):
        self.qc = qc
        self.sparse = sparse
        self.pp = pp
        self.genes = genes
        self.cells = cells
        self.annotate = annotate
        self.model = model
        self.pca = pca
        self.neighbors = neighbors
        self.clustering = clustering
        self.embedding = embedding

        if not self.clustering == None:
            self.neighbors = True
        if not self.neighbors:
            self.pca = True


    def __str__(self):
        return json.dumps(self.__dict__)
