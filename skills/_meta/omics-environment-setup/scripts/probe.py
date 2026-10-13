import importlib
import importlib.metadata as m
import json
import platform
packages = {'scanpy':'scanpy','anndata':'anndata','harmonypy':'harmonypy','python-igraph':'igraph','leidenalg':'leidenalg','umap-learn':'umap','numpy':'numpy','pandas':'pandas','scipy':'scipy','h5py':'h5py','pyyaml':'yaml','pillow':'PIL','pyarrow':'pyarrow','matplotlib-base':'matplotlib','scikit-misc':'skmisc','tifffile':'tifffile'}
distributions = {'python-igraph':'igraph','matplotlib-base':'matplotlib'}
versions = {'python':platform.python_version()}
for name, module in packages.items():
    importlib.import_module(module)
    versions[name] = m.version(distributions.get(name, name))
print(json.dumps({'versions': versions, 'machine':platform.machine()}))
