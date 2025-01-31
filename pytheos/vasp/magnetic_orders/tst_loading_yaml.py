import yaml

with open(f"rocksalt_prim_afmG.yaml", "r") as f:
    magmoms = yaml.load(f, Loader=yaml.SafeLoader)

print(magmoms)
