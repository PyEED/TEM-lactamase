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

### SetUp with Neo4J

This is the docker container in which data will be saved:
```
docker run -it --name pyeed-neo4j-niklas-tem \
  --user="$(id -u):$(id -g)" \
  -e NEO4J_AUTH=neo4j/12345678900 \
  -p 7474:7474 \
  -p 7687:7687 \
  -e NEO4J_PLUGINS='["graph-data-science", "apoc"]' \
  -e NEO4J_dbms_security_procedures_unrestricted="gds.*,apoc.*,gds.util.*" \
  -d neo4j:latest
```

