\# atelier-tensor-network



2025年度『ぶつりのあとりえ。』の記事「テンソルネットワークって？ by 沖縄のフォン・ノイマン」の図の作成に用いたプログラム．

プログラムはすべてJuliaで書かれている．



\# ファイルの説明と使い方



├───ctmrg

│       CTMRG\_main.jl    # CTMRG\[1]\[2]のアルゴリズム

│       ising\_ctmrg\_spm.jl    # 自発磁化のグラフの作成

│       ising\_ctmrg\_sscor.jl    # 最近接スピン相関のグラフの作成

│       ising\_dual\_ctmrg.jl    # 裏格子変換を行った場合の最近接スピン相関のグラフの作成

│

└───snapshot

&nbsp;       snapshot.jl    # スナップショット\[3]のアルゴリズム

&nbsp;       snapshot\_oil\_water.jl    # 水と油の図\[4]の作成

&nbsp;       snapshot\_something.jl    # テキトーに考えた図の作成



/ctmrg/CTMRG\_main.jl, /snapshot/snapshot.jlはアルゴリズムが書かれており，

実行しても何も起きない．

図の作成にはそれ以外の.jlファイル（snapshot/snapshot\_something.jlなど）を実行する．



\# 参考文献



\[1] T. Nishino and K. Okunishi, “Corner Transfer Matrix Renormalization Group

Method,” J. Phys. Soc. Jpn. 65, 891 (1996), arXiv:cond-mat/9507087.

\[2] T. Nishino and K. Okunishi, “Corner Transfer Matrix Algorithm for Classical Renor-

malization Group,” J. Phys. Soc. Jpn. 66, 3040 (1997), arXiv:cond-mat/9705072.

\[3] K. Ueda et al., “Snapshot Observation for 2D Classical Lattice Models by Cor-

ner Transfer Matrix Renormalization Group,” J. Phys. Soc. Jpn. 74, 111 (2005),

arXiv:cond-mat/0409445.

\[4]K. Ueda et al., “Snapshot Observation for 2D Classical Lattice Models by Cor-

ner Transfer Matrix Renormalization Group,” J. Phys. Soc. Jpn. 74, 111 (2005),

arXiv:cond-mat/0409445.

