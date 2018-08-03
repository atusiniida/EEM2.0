EEMを発現データA.tab、遺伝子セットデータA.gmtに適用
perl EEM2.0/perl/eem.pl  A.tab B.gmt 

A＿B.eemができる

各行にmodule の情報がはいっている

MICROTUBULE_CYTOSKELETON	[tab]3.820750485652578[tab]3[tab]	CoherenceBasedEEM absoluteRadius=0.34375 relativeRadius=0.050000000000000003 itrForPvalue2Calculation=300 coreGeneSize=10 seedGenes=[BIRC5, CDC2, DNALI1, FBXO5, LCK, MAP4, MAPT, NEK2, NUMA1, TOP2A] moduleGenes=[CDC2, NEK2, BIRC5, TOP2A] center=[-0.8829026374586828, 1.8518210320622268, 1.6803764005845103, 1.1561892297170608,.... 0.25079749074969349] Pvalue=4.2978498658743973 Pvalue1=3.0305692910087965 Pvalue2=4.2978498658743973

重要なところ
module ID
p-value  　
seedGenes: inputのgene set
moduleGenes: seedGenesのうち共発現しているもの


module clustering 
perl EEM2.0/perl/postEEMmoduleClustering.pl   A＿B.eem A.tab 

-log10 p値　X以上のモジュールをcorrelation cutoff Yでクラスタリング、module個数Z未満のクラスターは捨てる (Default: X = 6 Y=0.7 Z =0)
perl EEM2.0/perl/postEEMmoduleClustering.pl  -p X -c Y  -m Z A＿B.eem A.tab 



A＿B.posteem (ディレクトリ）ができる

A＿B.posteemのなかに
module.tab:  module activity profileのmatrix
module.pdf: module.tabのヒートマップ colored side barsがmodule clusterをしめす。
module.collapsed.tab:  module activity profileのクラスターをまとめたmeta module activity profile
module.collapsed.pdf: module.collapsed.tabのヒートマップ
module.gmt：module cluster中のmolule
gene.gmt: module cluster中の遺伝子 (module cluster中のmoduleのmodule geneをまとめたもの）
