# realtime-targeting-single-cell-long-read
Real-time targeted enrichment method for single-cell long-read sequencing.

### 📁 Repository Structure
```
realtime-targeting-scLongRead/
├── README.md                       # This file
├── envs/                           # Conda environments
├── resources/                      # Reference files and metadata
│   ├── annotations/                # GENCODE GTF and exon count files
│   ├── barcode2celltype/           # Barcode to cell type mapping 
│   ├── target_gene_sets/           # Target gene sets
│   ├── SRA-metadata.xlsx           # SRA metadata
│   └── celltype_mappings/          # Barcode to cell type
├── scripts/                        # Executable scripts
│   ├── count-umis-from-all-info.py # For UMI correction
│   ├── general/                    # General scripts
│   ├── scisorseq_scripts/          # Data preprocessing
│   ├── scisorATAC_scripts/         # Data preprocessing
│   └── isoquant_scripts/           # Data preprocessing
├── Realtime_targeting_workflow.md  # Targeting workflow
├── Analysis_main_figures.md        # Main analysis
└── Analysis_supplementary.md       # Supplementary analysis
```

### 💾 Data Availability
The raw sequencing data generated in this study have been deposited in the Sequence Read Archive database. 

