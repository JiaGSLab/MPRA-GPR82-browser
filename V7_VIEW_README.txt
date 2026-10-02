Claude v7 专用浏览包：hg38 / 170 bp insert / 2000 oligos

IGV: 解压后 File > Open Session 打开 igv/GPR82_v7_local.xml。请保留 tracks/ 相对目录。
在线会话: https://raw.githubusercontent.com/JiaGSLab/MPRA-GPR82-browser/main/igv/GPR82_v7_remote.xml
UCSC 主区: https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hubUrl=https%3A%2F%2Fraw.githubusercontent.com%2FJiaGSLab%2FMPRA-GPR82-browser%2Fmain%2Fhub_v7%2Fhub.txt&position=chrX%3A41703000-41751000
UCSC 全景: https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hubUrl=https%3A%2F%2Fraw.githubusercontent.com%2FJiaGSLab%2FMPRA-GPR82-browser%2Fmain%2Fhub_v7%2Fhub.txt&position=chrX%3A41435000-41755000
Hub URL: https://raw.githubusercontent.com/JiaGSLab/MPRA-GPR82-browser/main/hub_v7/hub.txt

自动带入项目选定的 cCRE、5 份巨噬细胞 DNase、巨噬细胞 H3K27ac、单核细胞 H3K27ac 代理、5 个条件 STAT3、GPR82 基因模型与 rE2G 预测连接。
古人类轨道是 Altai/Vindija/Denisova/Chagyrskaya 中通过本项目筛选的 26 个 ALT 位点注释，非原始古 DNA reads/BAM，也不是全基因组古人类轨道。
ENCODE 证据为项目缓存的局部峰/注释，不是自动加载 ENCODE 全部数据，不会随数据库更新自动刷新。峰条带不是测序信号强度曲线。
默认隐藏密集 v7 oligo 窗口和基因组对照；可在 UCSC 调整显示。IGV 可另加载 tracks/claude_v7_windows.bed 和 tracks/claude_v7_controls.bed。
原 v7 的 56 条 ladder 控制区间长度错误在本浏览包中按来源修正到 170 bp。原合成文件没有改动。
首次访问 UCSC 可能需要手动完成人机验证。
