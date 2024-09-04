# class or struct with list of processing options
# qc: bool
# sparce: bool
# pp: str (default, 'seurat', 'zheng17')
# pca: bool - required for neighbors
# neighbors: bool - required for clustering
# clustering: ('louvain', 'leiden')
# embedding: (['umap', 'tsne'])

import json

class Options:
    def __init__(
            self, 
            qc: bool = True,
            sparse: bool = True,
            # pp: str = 'default' | 'seurat' | 'zheng17',
            pp: str = 'default',
            pca: bool = True,
            neighbors: bool = True,
            # clustering: str = 'louvain'|'leiden',
            clustering: str = 'louvain',
            # embedding: list = ['umap', 'tsne']):
            embedding: list = ['umap']):
        self.qc = qc
        self.sparse = sparse
        self.pp = pp
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
