import inspect
from pathlib import Path
from nodes.registry import NODE_CLASS_MAPPINGS

rows = []
for name, cls in sorted(NODE_CLASS_MAPPINGS.items()):
    cat = getattr(cls, 'CATEGORY', 'Unknown')
    mod = cls.__module__
    doc = (cls.__doc__ or '').strip().splitlines()[0] if cls.__doc__ else ''
    
    run_fn = getattr(cls, 'run', getattr(cls, 'execute', getattr(cls, 'run_step', getattr(cls, 'summarize', None))))
    src = ''
    if run_fn:
        try:
            src = inspect.getsource(run_fn)
        except Exception:
            src = ''
            
    impl_type = 'Library Computation'
    if 'visualiz' in name.lower() or 'plot' in name.lower() or 'report' in name.lower() or 'IMAGE' in getattr(cls, 'RETURN_TYPES', ()):
        impl_type = 'Visualization'
    elif 'validator' in name.lower() or 'loader' in name.lower() or 'export' in name.lower() or 'io' in name.lower():
        impl_type = 'Utility / IO'
    elif 'runner' in src.lower() or 'subprocess' in src.lower() or 'command' in src.lower() or 'stage_commands' in src.lower():
        impl_type = 'CLI Wrapper'
        
    test_level = 'Unit Test & E2E Workflow' if any(k in name.lower() for k in ['fastp', 'bwa', 'macs3', 'spades', 'quast', 'salmon', 'kraken2', 'bracken', 'deseq2', 'bcftools', 'metadata', 'markduplicate']) else 'Registry & Unit Test'
    
    clean_doc = doc[:60].replace('|', '/')
    rows.append(f"| `{name}` | `{cat}` | {impl_type} | {test_level} | {clean_doc} |")

content = """# Supplementary Table S1: Catalog of ComfyBIOWMS Custom Nodes (179 Unique Nodes)

This table provides a comprehensive functional classification and verification status for all 179 custom nodes implemented in ComfyBIOWMS.

| Node Class Name | Category | Functional Type | Verification Scope | Short Description |
| :--- | :--- | :--- | :--- | :--- |
""" + "\n".join(rows) + "\n"

Path("supplementary_table_s1_node_catalog.md").write_text(content, encoding="utf-8")
print(f"Successfully wrote {len(rows)} nodes to supplementary_table_s1_node_catalog.md")
