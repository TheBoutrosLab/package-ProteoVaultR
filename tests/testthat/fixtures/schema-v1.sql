CREATE TABLE metadata (key VARCHAR PRIMARY KEY, value VARCHAR NOT NULL);
CREATE TABLE analysis (
    analysis_id BIGINT PRIMARY KEY, dataset_id VARCHAR NOT NULL,
    pipeline VARCHAR NOT NULL, pipeline_version VARCHAR, root_path VARCHAR NOT NULL
);
CREATE TABLE source_files (
    source_id BIGINT PRIMARY KEY, path VARCHAR NOT NULL, relative_path VARCHAR NOT NULL,
    format VARCHAR NOT NULL, stage VARCHAR NOT NULL, size_bytes BIGINT NOT NULL,
    mtime_ns BIGINT NOT NULL, status VARCHAR NOT NULL, details VARCHAR
);
CREATE TABLE samples (sample_id BIGINT PRIMARY KEY, sample_name VARCHAR NOT NULL);
CREATE TABLE runs (
    run_id BIGINT PRIMARY KEY, sample_id BIGINT,
    run_name VARCHAR NOT NULL, file_name VARCHAR
);
CREATE TABLE spectra (
    spectrum_id BIGINT PRIMARY KEY, run_id BIGINT NOT NULL,
    native_id VARCHAR NOT NULL, rt_seconds DOUBLE, ms_level INTEGER,
    precursor_mz DOUBLE, precursor_charge INTEGER, source_id BIGINT,
    CHECK (length(native_id) > 0)
);
CREATE TABLE peptides (peptide_id BIGINT PRIMARY KEY, sequence VARCHAR NOT NULL);
CREATE TABLE peptidoforms (
    peptidoform_id BIGINT PRIMARY KEY, peptide_id BIGINT NOT NULL,
    modified_sequence VARCHAR NOT NULL, notation VARCHAR NOT NULL,
    CHECK (length(modified_sequence) > 0)
);
CREATE TABLE precursors (
    precursor_id BIGINT PRIMARY KEY,
    peptidoform_id BIGINT NOT NULL, charge INTEGER NOT NULL CHECK(charge > 0)
);
CREATE TABLE psms (
    psm_id BIGINT PRIMARY KEY, spectrum_id BIGINT NOT NULL,
    precursor_id BIGINT NOT NULL, search_engine VARCHAR NOT NULL,
    candidate_rank INTEGER NOT NULL, target_decoy VARCHAR, source_id BIGINT NOT NULL,
    identification_index INTEGER NOT NULL, hit_index INTEGER NOT NULL, stage_id INTEGER NOT NULL
);
CREATE TABLE psm_scores (
    psm_id BIGINT NOT NULL, score_id INTEGER NOT NULL, value DOUBLE NOT NULL,
    higher_better BOOLEAN, PRIMARY KEY(psm_id, score_id)
);
CREATE TABLE psm_features (
    psm_id BIGINT NOT NULL, name VARCHAR NOT NULL, value DOUBLE NOT NULL,
    PRIMARY KEY(psm_id, name)
);
CREATE TABLE psm_metadata (
    psm_id BIGINT NOT NULL, name VARCHAR NOT NULL, value VARCHAR NOT NULL,
    PRIMARY KEY(psm_id, name)
);
CREATE TABLE proteins (
    protein_id BIGINT PRIMARY KEY, accession VARCHAR NOT NULL,
    sequence VARCHAR, description VARCHAR, is_decoy BOOLEAN, is_contaminant BOOLEAN
);
CREATE TABLE protein_identifiers (
    protein_id BIGINT NOT NULL, source_id BIGINT NOT NULL,
    source_accession VARCHAR NOT NULL, mapping_type VARCHAR NOT NULL,
    mapped_id VARCHAR NOT NULL,
    PRIMARY KEY (protein_id, source_id, mapping_type, mapped_id)
);
CREATE TABLE peptide_proteins (
    peptide_id BIGINT NOT NULL, protein_id BIGINT NOT NULL,
    PRIMARY KEY(peptide_id, protein_id)
);
CREATE TABLE psm_proteins (
    psm_id BIGINT NOT NULL, protein_id BIGINT NOT NULL,
    PRIMARY KEY(psm_id, protein_id)
);
CREATE TABLE protein_groups (
    group_id BIGINT PRIMARY KEY, source_id BIGINT NOT NULL,
    source_group_id VARCHAR NOT NULL, method_id INTEGER NOT NULL, group_type VARCHAR NOT NULL,
    score DOUBLE, score_type VARCHAR
);
CREATE TABLE protein_group_members (
    group_id BIGINT NOT NULL, protein_id BIGINT NOT NULL,
    PRIMARY KEY(group_id, protein_id)
);
CREATE TABLE quant_channels (
    channel_id BIGINT PRIMARY KEY, sample_id BIGINT NOT NULL,
    run_id BIGINT, label VARCHAR NOT NULL
);
CREATE TABLE precursor_quant (
    precursor_id BIGINT NOT NULL, channel_id BIGINT NOT NULL,
    source_id BIGINT NOT NULL, method_id INTEGER NOT NULL,
    quantity DOUBLE, q_value DOUBLE, pep DOUBLE,
    PRIMARY KEY(precursor_id, channel_id, source_id, method_id)
);
CREATE TABLE protein_quant (
    group_id BIGINT NOT NULL, channel_id BIGINT NOT NULL,
    source_id BIGINT NOT NULL, method_id INTEGER NOT NULL,
    quantity DOUBLE, PRIMARY KEY(group_id, channel_id, source_id, method_id)
);
CREATE TABLE precursor_groups (
    precursor_id BIGINT NOT NULL, group_id BIGINT NOT NULL,
    PRIMARY KEY(precursor_id, group_id)
);
CREATE TABLE pin_features (
    source_id BIGINT NOT NULL, row_number INTEGER NOT NULL,
    spectrum_reference VARCHAR NOT NULL, scan_number VARCHAR, label INTEGER,
    feature_name VARCHAR NOT NULL, value DOUBLE NOT NULL,
    PRIMARY KEY(source_id, row_number, feature_name)
);
CREATE TABLE qc_metrics (
    source_id BIGINT NOT NULL,
    row_number BIGINT NOT NULL,
    dataset_id VARCHAR NOT NULL,
    layer VARCHAR NOT NULL,
    source_process VARCHAR,
    sample_name VARCHAR,
    search_engine VARCHAR,
    fraction_id VARCHAR,
    metric_name VARCHAR NOT NULL,
    "group" VARCHAR,
    value VARCHAR,
    unit VARCHAR,
    status VARCHAR NOT NULL,
    source_file VARCHAR,
    PRIMARY KEY(source_id, row_number)
);
CREATE TABLE import_counts (
    source_id BIGINT NOT NULL, entity VARCHAR NOT NULL,
    row_count BIGINT NOT NULL, PRIMARY KEY(source_id, entity)
);
CREATE TABLE quant_methods (method_id INTEGER PRIMARY KEY, name VARCHAR NOT NULL);
CREATE TABLE grouping_methods (method_id INTEGER PRIMARY KEY, name VARCHAR NOT NULL);
CREATE TABLE processing_stages (stage_id INTEGER PRIMARY KEY, name VARCHAR NOT NULL);
CREATE TABLE score_types (score_id INTEGER PRIMARY KEY, name VARCHAR NOT NULL);
CREATE TABLE score_aliases (score_id INTEGER NOT NULL, source_name VARCHAR NOT NULL, PRIMARY KEY(score_id, source_name));
CREATE VIEW protein_groups_readable AS SELECT g.* EXCLUDE(method_id), t.name AS method FROM protein_groups g JOIN grouping_methods t USING(method_id);
CREATE VIEW psm_scores_readable AS SELECT q.psm_id, t.name, q.value, q.higher_better FROM psm_scores q JOIN score_types t USING(score_id);
CREATE VIEW psms_readable AS SELECT p.* EXCLUDE(stage_id), t.name AS stage FROM psms p JOIN processing_stages t USING(stage_id);
CREATE VIEW protein_quant_readable AS SELECT q.* EXCLUDE(method_id), t.name AS method FROM protein_quant q JOIN quant_methods t USING(method_id);
CREATE VIEW precursor_quant_readable AS SELECT q.* EXCLUDE(method_id), t.name AS method FROM precursor_quant q JOIN quant_methods t USING(method_id);
CREATE VIEW protein_quant_long AS
    SELECT s.sample_name, g.group_id, g.source_group_id, gm.name AS grouping_method,
        t.name AS quantification_method, q.quantity, q.source_id
    FROM protein_quant q JOIN quant_methods t USING(method_id) JOIN protein_groups g USING(group_id) JOIN grouping_methods gm ON g.method_id=gm.method_id
    JOIN quant_channels c USING(channel_id) JOIN samples s USING(sample_id);
CREATE VIEW precursor_quant_long AS
    SELECT s.sample_name, p.precursor_id, f.modified_sequence, f.notation, p.charge,
        t.name AS method, q.quantity, q.q_value, q.pep, q.source_id
    FROM precursor_quant q JOIN quant_methods t USING(method_id) JOIN precursors p USING(precursor_id)
    JOIN peptidoforms f USING(peptidoform_id) JOIN quant_channels c USING(channel_id)
    JOIN samples s USING(sample_id);
CREATE VIEW protein_group_summary AS
    SELECT g.*, count(m.protein_id) AS member_count FROM protein_groups_readable g
    LEFT JOIN protein_group_members m USING(group_id) GROUP BY ALL;
CREATE VIEW sample_protein_counts AS
    SELECT sample_name, grouping_method, quantification_method, source_id,
        count(*) FILTER (WHERE quantity > 0) AS quantified_protein_groups
    FROM protein_quant_long GROUP BY ALL;
