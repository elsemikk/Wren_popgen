I will view the Ficedula genome sequence here: https://www.ncbi.nlm.nih.gov/genome/gdv/browser/genome/?id=GCF_000247815.1
I am interested in the order of genes in the segment at 17.7-22.4 Mb. Actually, I should be looking 2 Mb up and down stream to check synteny there!
 so under "search assembly", I searched `8:15,700,000-24,400,000`.
Then, I obtained the data as a Genbank flat file with Download -> Genbank Flat File -> Visible range

```bash
cd ~/Documents/GitHub/Wren_popgen/
#get a list of genes in this region
cat NC_021680.1[15700000..24400000].flat | grep "/gene" | sed 's/\/gene=//g' | sed 's/"//g' | sed 's/^ *//g' | grep -v "gene_synonym" > Ficedula_chr8genes
#get a list of their locations
cat NC_021680.1[15700000..24400000].flat | grep "  gene  " | sed 's/^ *gene *//g' > Ficedula_chr8loci
```

Next I looked up Taeniopygia here: https://www.ncbi.nlm.nih.gov/genome/gdv/?org=taeniopygia-guttata&group=passeriformes
I searched for the first gene in the Ficedula segment, which was HS2ST1. Then I looked at this 4.7 Mb region starting at the start of that gene, plus 2 Mb flanking, which would be `8:10,700,000-23,400,000`. Then I downloaded this data again.
```bash
#get a list of genes in this region
cat NC_045007.1[12700000..21400000].flat | grep "/gene" | sed 's/\/gene=//g' | sed 's/"//g' | sed 's/^ *//g' | grep -v "gene_synonym" > Taeniopygia_chr8genes
#get a list of their locations
cat NC_045007.1[12700000..21400000].flat | grep "  gene  " | sed 's/^ *gene *//g' > Taeniopygia.chr8loci
```

Lepidothrix coronata
HS2ST1 is at 9.0 Mb, so get 7-15.7 Mb, however the scaffold is too short, only got NW_016690192.1:7,000,000 - 10,516,705
```bash

#get a list of genes in this region
cat NW_016690192.1[7000000..10516705].flat | grep "/gene" | sed 's/\/gene=//g' | sed 's/"//g' | sed 's/^ *//g' | grep -v "gene_synonym" > Lepidothrix.chr8genes
#get a list of their locations
cat NW_016690192.1[7000000..10516705].flat | grep "  gene  " | sed 's/^ *gene *//g' > Lepidothrix.chr8loci
http://127.0.0.1:39419/graphics/plot_zoom_png?width=1200&height=835
#get a list of genes in this region
cat NW_016690316.1[1..2590296].flat | grep "/gene" | sed 's/\/gene=//g' | sed 's/"//g' | sed 's/^ *//g' | grep -v "gene_synonym" > Lepidothrix2.chr8genes
#get a list of their locations
cat NW_016690316.1[1..2590296].flat | grep "  gene  " | sed 's/^ *gene *//g' > Lepidothrix2.chr8loci

```

I realized that Manacus is in reverse complementary sequence compared to Ficedula. HS2ST1 is at 2.4 Mb. So the sequence will be a fragment.

Manacus vitellinus: `NW_021939546.1:1-4,400,000`
```bash
#get a list of genes in this region
cat NW_021939546.1[400000..7100000].flat | grep "/gene" | sed 's/\/gene=//g' | sed 's/"//g' | sed 's/^ *//g' | grep -v "gene_synonym" > Manacus.chr8genes
#get a list of their locations
cat NW_021939546.1[400000..7100000].flat | grep "  gene  " | sed 's/^ *gene *//g' > Manacus.chr8loci

#get a list of genes in this region
cat NW_021940660.1[1..3820261].flat | grep "/gene" | sed 's/\/gene=//g' | sed 's/"//g' | sed 's/^ *//g' | grep -v "gene_synonym" > Manacus2.chr8genes
#get a list of their locations
cat NW_021940660.1[1..3820261].flat | grep "  gene  " | sed 's/^ *gene *//g' > Manacus2.chr8loci

```

Empidonax
```bash
#get a list of genes in this region
cat NW_020955307.1[1..1751862].flat | grep "/gene" | sed 's/\/gene=//g' | sed 's/"//g' | sed 's/^ *//g' | grep -v "gene_synonym" > Empidonax.chr8genes
#get a list of their locations
cat NW_020955307.1[1..1751862].flat | grep "  gene  " | sed 's/^ *gene *//g' > Empidonax.chr8loci
```

Gallus: 16.2, reverse -> 8:1,200,000-17,700,000

```bash
#get a list of genes in this region
cat NC_006095.5[1200000..17700000].flat | grep "/gene" | sed 's/\/gene=//g' | sed 's/"//g' | sed 's/^ *//g' | grep -v "gene_synonym" > Gallus.chr8genes
#get a list of their locations
cat NC_006095.5[1200000..17700000].flat | grep "  gene  " | sed 's/^ *gene *//g' > Gallus.chr8loci
```
Corvus
```bash
#get a list of genes in this region
cat NW_008237239.1[4900000..13100000].flat | grep "/gene" | sed 's/\/gene=//g' | sed 's/"//g' | sed 's/^ *//g' | grep -v "gene_synonym" > Corvus.chr8genes
#get a list of their locations
cat NW_008237239.1[4900000..13100000].flat | grep "  gene  " | sed 's/^ *gene *//g' > Corvus.chr8loci
```

Discarded:
Lonchura striata: `8:14,200,000-20,900,000`
```bash
#get a list of genes in this region
cat NC_042574.1[14200000..20900000].flat | grep "/gene" | sed 's/\/gene=//g' | sed 's/"//g' | sed 's/^ *//g' | grep -v "gene_synonym" > Lonchura.chr8genes
#get a list of their locations
cat NC_042574.1[14200000..20900000].flat | grep "  gene  " | sed 's/^ *gene *//g' > Lonchura.chr8loci
```
