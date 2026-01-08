#defining directory with DMR bed files
dmrdir='/DMR/BEDs/'

#defining output directory
outdir='/output/directory/'

#running Homer to assess enrichment for transcription factor binding sites among top 2000 hypermethylated and hypomethylated region
Homer findMotifsGenome.pl $dmrdir/Female.hyper2000.bed hg38 $outdir -size 300 -useNewBg -p 10

Homer findMotifsGenome.pl $dmrdir/Female.hypo2000.bed hg38 $outdir -size 300 -useNewBg -p 10

Homer findMotifsGenome.pl $dmrdir/Male.hyper2000.bed hg38 $outdir -size 300 -useNewBg -p 10

Homer findMotifsGenome.pl $dmrdir/Male.hypo2000.bed hg38 $outdir -size 300 -useNewBg -p 10

