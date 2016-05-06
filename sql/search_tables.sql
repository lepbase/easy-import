DROP TABLE IF EXISTS multi;
CREATE TABLE IF NOT EXISTS multi (
  feature_id               INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
  feature_type             ENUM('gene','seq_region','species') NOT NULL,
  core_id                  INT(10) UNSIGNED NOT NULL,
  string                   VARCHAR(255) NOT NULL,
  string_type              VARCHAR(128) NOT NULL,
  detail              	   VARCHAR(255),
  production_name          VARCHAR(128) NOT NULL NOT NULL,
   
  PRIMARY KEY (feature_id),
  KEY core_idx (core_id),
  KEY pn_feat_core_idx (production_name, feature_type, core_id),
  KEY str_idx (string),
  KEY pn_str_idx (production_name, string),
  KEY pn_feat_str_idx (production_name, feature_type, string)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;

DROP TABLE IF EXISTS multi_255;
CREATE TABLE IF NOT EXISTS multi_255 (
  feature_id               INT(10) UNSIGNED NOT NULL,
  feature_type             ENUM('gene','seq_region','species') NOT NULL,
  core_id                  INT(10) UNSIGNED NOT NULL,
  string                   VARCHAR(255) NOT NULL,
  string_type              VARCHAR(128) NOT NULL,
  production_name          VARCHAR(128) NOT NULL NOT NULL,
   
  PRIMARY KEY (feature_id),
  KEY core_idx (core_id),
  KEY pn_feat_core_idx (production_name, feature_type, core_id),
  KEY str_idx (string),
  KEY pn_str_idx (production_name, string),
  KEY pn_feat_str_idx (production_name, feature_type, string)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;
ALTER TABLE multi_255 ADD fulltext index string_index (string);

DROP TABLE IF EXISTS multi_32;
CREATE TABLE IF NOT EXISTS multi_32 (
  feature_id               INT(10) UNSIGNED NOT NULL,
  feature_type             ENUM('gene','seq_region','species') NOT NULL,
  core_id                  INT(10) UNSIGNED NOT NULL,
  string                   VARCHAR(32) NOT NULL,
  string_type              VARCHAR(128) NOT NULL,
  production_name          VARCHAR(128) NOT NULL NOT NULL,
   
  PRIMARY KEY (feature_id),
  KEY core_idx (core_id),
  KEY pn_feat_core_idx (production_name, feature_type, core_id),
  KEY str_idx (string),
  KEY pn_str_idx (production_name, string),
  KEY pn_feat_str_idx (production_name, feature_type, string)
) COLLATE=latin1_swedish_ci ENGINE=MyISAM;


