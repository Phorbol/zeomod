from collections import Counter
from itertools import combinations
import numpy as np
from ase import Atoms

def validate_pure_silica(atoms: Atoms):
    symbols = atoms.get_chemical_symbols()
    counts = Counter(symbols)
    # 电荷守恒: 4*Si + 1*H - 2*O = 0
    charge = 4*counts['Si'] + 1*counts['H'] - 2*counts['O']
    if charge == 0:
        print(f"✅ 纯硅骨架验证通过 (Si{counts['Si']}O{counts['O']}H{counts['H']})")
    else:
        print(f"❌ 纯硅骨架电荷不平衡! 净电荷: {charge}")

def validate_bronsted(atoms: Atoms, graph):
    si = len([a for a in atoms if a.symbol == 'Si'])
    al = len([a for a in atoms if a.symbol == 'Al'])
    o = len([a for a in atoms if a.symbol == 'O'])
    h = len([a for a in atoms if a.symbol == 'H'])
    
    # 1. 电荷
    charge = 4*si + 3*al + 1*h - 2*o
    if charge != 0:
        print(f"❌ 最终结构电荷不平衡! Net: {charge}")
    else:
        print("✅ 电荷平衡验证通过。")
        
    # 2. Löwenstein
    al_indices = [a.index for a in atoms if a.symbol == 'Al']
    errors = 0
    for i, j in combinations(al_indices, 2):
        if graph.get_distance(i, j) <= 1:
            print(f"❌ 违反 Löwenstein 规则: Al({i})-Al({j})")
            errors += 1
    
    if errors == 0:
        print("✅ Löwenstein 规则验证通过 (无 Al-O-Al)。")