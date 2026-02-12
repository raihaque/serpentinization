# serpentinization
modeling olivine to serpentinization using lmgc90

-> create_hex_groups.py:
    Creates groups of hexagons
    Stores hexagon data in a json file

-> gen_olivine.py:
    Reads hexagon data from the json file
    Generates the geometry for pylmgc90,

-> command.py:
    Runs main simulation
