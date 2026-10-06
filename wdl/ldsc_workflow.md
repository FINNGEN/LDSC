# ldsc_rg Workflow

```mermaid
flowchart TD
    IN["INPUTS\nmeta_fg · meta_other · population · name\nonly_het · snplist · ld_root"]
    FM["filter_meta\nsplit both meta files into chunks"]
    D1{"meta_fg\n≠ meta_other?"}

    subgraph SFG["↻  scatter — fg chunks"]
        PFG["premunge_fg\nRSID map · column standardise"]
        MFG["munge_ldsc  fg\nmunge_sumstats · heritability"]
        PFG --> MFG
    end

    subgraph SOT["↻  scatter — other chunks"]
        POT["premunge_other\nRSID map · column standardise"]
        MOT["munge_ldsc  other\nmunge_sumstats · heritability"]
        POT --> MOT
    end

    COMB["flatten / combine\nall_het_jsons · all_het_logs · all_munged"]
    H2["gather_h2\nextract_metadata · plot_summary"]
    D2{"only_het\n= false?"}
    RC["return_couples\nall phenotype pairs → chunked couples"]

    subgraph SMRG["↻  scatter — couple chunks"]
        MRG["multi_rg\nldsc_mult.py · genetic correlations"]
    end

    GS["gather_summaries\nmerge logs · concat summary TSVs"]
    OUT["OUTPUTS\nherit_tsv · herit_log · munged_ss[]\ncorr_summary · corr_log"]

    IN --> FM
    FM -->|"chunk_fg[]"| SFG
    FM -->|"chunk_other[]"| D1
    D1 -->|"Yes"| SOT
    D1 -->|"No"| COMB
    MFG -->|"het_json · het_log · munged[]"| COMB
    MOT -->|"het_json · het_log · munged[]"| COMB
    COMB -->|"het_jsons · het_logs"| H2
    COMB -->|"munged_ss[]"| OUT
    COMB --> D2
    D2 -->|"Yes  munged[]"| RC
    D2 -->|"No"| OUT
    RC -->|"couples[] · paths_list[]"| SMRG
    MRG -->|"summary[] · log[]"| GS
    H2 -->|"herit_tsv · herit_log"| OUT
    GS -->|"corr_summary · corr_log"| OUT

    style IN   fill:#4a90d9,color:#fff,stroke:#3a80c9
    style OUT  fill:#4a90d9,color:#fff,stroke:#3a80c9
    style FM   fill:#2ecc71,color:#fff,stroke:#27ae60
    style PFG  fill:#2ecc71,color:#fff,stroke:#27ae60
    style MFG  fill:#2ecc71,color:#fff,stroke:#27ae60
    style COMB fill:#7f8c8d,color:#fff,stroke:#636e72
    style H2   fill:#2ecc71,color:#fff,stroke:#27ae60
    style POT  fill:#9b59b6,color:#fff,stroke:#8e44ad
    style MOT  fill:#9b59b6,color:#fff,stroke:#8e44ad
    style RC   fill:#e67e22,color:#fff,stroke:#d35400
    style MRG  fill:#e67e22,color:#fff,stroke:#d35400
    style GS   fill:#e67e22,color:#fff,stroke:#d35400
    style D1   fill:#ecf0f1,stroke:#95a5a6
    style D2   fill:#ecf0f1,stroke:#95a5a6
```

## Colour key

| Colour | Meaning |
|--------|---------|
| Blue | Workflow inputs / outputs |
| Green | Unconditional task |
| Purple | Conditional task — runs only when `meta_fg ≠ meta_other` |
| Orange | Conditional task — runs only when `only_het = false` |
| Gray | Flatten / combine step |
```
