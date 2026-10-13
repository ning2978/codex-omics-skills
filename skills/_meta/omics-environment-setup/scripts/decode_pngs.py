from pathlib import Path
import sys
from PIL import Image
paths = list(Path(sys.argv[1]).rglob('*.png'))
assert paths, 'No PNG outputs'
for path in paths:
    with Image.open(path) as image:
        image.load()
        assert image.width > 0 and image.height > 0
print(f'PASS decoded {len(paths)} PNGs')
