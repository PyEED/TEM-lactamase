# PyEED Data Model

This document describes the data model for PyEED, representing structured nodes and relationships for DNA, protein, ontology objects, and annotations. The model leverages a graph-based structure to support genotype-phenotype exploration and integration with biological ontologies.

---

### Organism

This object represents an organism associated with the DNA and protein sequences in the database.

- **taxonomy_id**
  - Type: integer
  - Description: Unique taxonomy ID for the organism.
  - PK: True
- **name**
  - Type: string
  - Description: Scientific name of the organism.

---

### DNA

- **accession_id**
  - Type: string
  - Description: Unique identifier for the DNA sequence.
  - PK: True
- **sequence**
  - Type: string
  - Description: The nucleotide sequence of the DNA.
- **name**
  - Type: string
  - Description: Name or description of the DNA sequence.
- **seq_length**
  - Type: integer
  - Description: Length of the DNA sequence.
- **gc_content**
  - Type: float
  - Description: GC content percentage of the DNA.
- **organism**
  - Type: [Organism](#organism)
  - Description: The organism to which this DNA belongs.
- **proteins**
  - Type: [Protein](#protein)[]
  - Description: Proteins encoded by this DNA.
- **sites**
  - Type: [Site](#site)[]
  - Description: Sites associated with this DNA.

---

### Protein

- **accession_id**
  - Type: string
  - Description: Unique identifier for the protein sequence.
  - PK: True
- **sequence**
  - Type: string
  - Description: The amino acid sequence of the protein.
- **name**
  - Type: string
  - Description: Name or description of the protein.
- **seq_length**
  - Type: integer
  - Description: Length of the protein sequence.
- **mol_weight**
  - Type: float
  - Description: Molecular weight of the protein.
- **ec_number**
  - Type: string
  - Description: Enzyme Commission number.
- **dna**
  - Type: [DNA](#dna)
  - Description: The DNA sequence that encodes this protein.
- **sites**
  - Type: [Site](#site)[]
  - Description: Sites associated with this protein.
- **regions**
  - Type: [Region](#region)[]
  - Description: Regions associated with this protein.
- **mutations**
  - Type: [Mutation](#mutation)[]
  - Description: Mutations involving this protein.

---

### Site

- **site_id**
  - Type: string
  - Description: Unique identifier for the site.
  - PK: True
- **name**
  - Type: string
  - Description: Name of the site.
- **annotation**
  - Type: string
  - Description: Type of site (e.g., binding site, active site).
- **positions**
  - Type: array[integer]
  - Description: Positions in the DNA or protein sequence associated with the site.
- **protein**
  - Type: [Protein](#protein)
  - Description: The protein where this site is found.

---

### Region

- **region_id**
  - Type: string
  - Description: Unique identifier for the region.
  - PK: True
- **annotation**
  - Type: string
  - Description: Annotation type of the region (e.g., domain, motif).
- **start**
  - Type: integer
  - Description: Start position of the region.
- **end**
  - Type: integer
  - Description: End position of the region.
- **protein**
  - Type: [Protein](#protein)
  - Description: The protein this region belongs to.

---

### StandardNumbering

- **name**
  - Type: string
  - Description: Name of the standard numbering schema.
- **definition**
  - Type: string
  - Description: Detailed description of the schema.
- **proteins**
  - Type: [Protein](#protein)[]
  - Description: Proteins aligned to this standard numbering schema.

---

### GOAnnotation

- **go_id**
  - Type: string
  - Description: Gene Ontology ID.
  - PK: True
- **term**
  - Type: string
  - Description: Name of the GO term.
- **definition**
  - Type: string
  - Description: Definition of the GO term.
- **proteins**
  - Type: [Protein](#protein)[]
  - Description: Proteins associated with this GO annotation.

---

### Mutation

- **from_positions**
  - Type: array[integer]
  - Description: Positions in the original sequence where mutations occur.
- **to_positions**
  - Type: array[integer]
  - Description: Positions in the mutated sequence corresponding to the mutation.
- **from_monomers**
  - Type: array[string]
  - Description: Original residues or nucleotides at the mutation positions.
- **to_monomers**
  - Type: array[string]
  - Description: Mutated residues or nucleotides at the mutation positions.
- **original_sequence**
  - Type: [Protein](#protein) | [DNA](#dna)
  - Description: The original sequence.
- **mutated_sequence**
  - Type: [Protein](#protein) | [DNA](#dna)
  - Description: The mutated sequence.

---

### OntologyObject

- **name**
  - Type: string
  - Description: Name of the ontology object.
  - PK: True
- **description**
  - Type: string
  - Description: Description of the ontology object.
- **label**
  - Type: string
  - Description: Label for the ontology object.
- **synonyms**
  - Type: array[string]
  - Description: Synonyms for the ontology object.
- **subclasses**
  - Type: [OntologyObject](#ontologyobject)[]
  - Description: Subclasses of this ontology object.
- **custom_relationships**
  - Type: array[CustomRelationship](#customrelationship)
  - Description: Custom relationships involving this ontology object.
