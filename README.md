# TEM-lactamase

#### TO-DO NEXT
- [ ] starte mit der Proteinsequenz von TEM-1
- [ ] BLAST -> 1000-2000  ähnliche Sequenzen
- [ ] Multiple sequence alignment (Clustal omega)
- [ ] Annotation in TEM-1: signal sequence, TEM domain, C-terminus (alles, was C-terminal an die "TEM domain" anschließt, z.B. His-tag)
- [ ] Übertragung der Annotation auf alle Homologe
- [ ] Ermittelung der TEM-number aus dem Proteinname: TEM-1, TEM-2,...
- [ ] Filtern der "TEM domain": identische Sequenzen entfernen; Auswahl: die Sequenz behalten, deren TEM-number bekannt ist
- [ ] Multiple pairwise alignment der "TEM domain", Ermittelung der Anzahl der Aminosäureaustausche zwischen den Sequenzen
- [ ] Netzwerk: edge = 1 Aminosäureaustausch, labeling mit TEM-number 


enroot import docker://haeussma/pyeed-notebook:latest
enroot create --name haeussma+pyeed-notebook+latest+new haeussma+pyeed-notebook+latest.sqsh
