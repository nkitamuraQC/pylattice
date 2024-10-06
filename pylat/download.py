from pymatgen.ext.matproj import MPRester

def download_via_pymatgen(material_id):
    API_KEY = "YOUR_API_KEY"
    
    # MPResterインスタンスを作成
    with MPRester(API_KEY) as m:
        # 取得したい材料のMaterials Project IDを指定
        material_id = material_id  
    
        # 構造データを取得
        structure = m.get_structure_by_material_id(material_id)
    
        # CIFファイルとして保存
        structure.to(fmt="cif", filename=f"{material_id}.cif")
    return
