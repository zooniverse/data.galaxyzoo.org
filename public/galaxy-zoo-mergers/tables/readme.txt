Columns for target_info.txt

ID	Target id, one-up number roughly same as Merger Zoo presentation order
SDSSID	SDSS DR7 ID for target
PRI_RA_DEG	RA in degrees for primary galaxy
PRI_DEC_DEG	DEC in degrees for primary galaxy
PRI_NAMES	Alternate names for primary galaxy: SDSS, Arp, NGC/IC, UGC, 2MASS
SEC_RA_DEG	RA in degrees for secondary galaxy
SEC_DEC_DEG	DEC in degrees for secondary galaxy
SEC_NAMES	Alternate names for secondary galaxy: SDSS, Arp, NGC/IC, UGC, 2MASS



Columns for target_sdss_dr7_dr8.txt
For DR8 values, the ra/dec from DR7 was used with a call to the CAS database function dbo.fGetNearbyObjEq() using a radius of 0.75 arcmin.
Target 54 was from DR8 originally, no corresponding values for DR7 found.  For several DR7 targets no DR8 values were found.

ID	Target id, one-up number roughly same as Merger Zoo presentation order
SDSSID	SDSS DR7 ID for target
PRI_DR7_ID	SDSS DR7 id for primary galaxy
PRI_DR8_ID	SDSS DR8 id for primary galaxy
SEC_DR7_ID	SDSS DR7 id for secondary galaxy
SEC_DR8_ID	SDSS DR8 id for secondary galaxy
PRI_DR7_U	Dered U magnitude (model - extinction) from DR7 for primary
PRI_DR7_G	Dered G magnitude (model - extinction) from DR7 for primary
PRI_DR7_R	Dered R magnitude (model - extinction) from DR7 for primary
PRI_DR7_I	Dered I magnitude (model - extinction) from DR7 for primary
PRI_DR7_Z	Dered Z magnitude (model - extinction) from DR7 for primary
PRI_DR7_SPECZ	Spectral redshift from DR7 for primary
PRI_DR7_PHOTOZ	Photometric redshift from DR7 for primary
PRI_DR7_PHOTOZ2	Alternate photometric redshift from DR7 for primary
PRI_DR8_U	Dered U magnitude (model - extinction) from DR8 for primary
PRI_DR8_G	Dered G magnitude (model - extinction) from DR8 for primary
PRI_DR8_R	Dered R magnitude (model - extinction) from DR8 for primary
PRI_DR8_I	Dered I magnitude (model - extinction) from DR8 for primary
PRI_DR8_Z	Dered Z magnitude (model - extinction) from DR8 for primary
PRI_DR8_SPECZ	Spectral redshift from DR8 for primary
PRI_DR8_PHOTOZ	Photometric redshift from DR8 for primary
PRI_DR8_PHOTOZ2	Alternate photometric redshift from DR8 for primary
SEC_DR7_U	Dered U magnitude (model - extinction) from DR7 for secondary
SEC_DR7_G	Dered G magnitude (model - extinction) from DR7 for secondary
SEC_DR7_R	Dered R magnitude (model - extinction) from DR7 for secondary
SEC_DR7_I	Dered I magnitude (model - extinction) from DR7 for secondary
SEC_DR7_Z	Dered Z magnitude (model - extinction) from DR7 for secondary
SEC_DR7_SPECZ	Spectral redshift from DR7 for secondary
SEC_DR7_PHOTOZ	Photometric redshift from DR7 for secondary
SEC_DR7_PHOTOZ2	Alternate photometric redshift from DR7 for secondary
SEC_DR8_U	Dered U magnitude (model - extinction) from DR8 for secondary
SEC_DR8_G	Dered G magnitude (model - extinction) from DR8 for secondary
SEC_DR8_R	Dered R magnitude (model - extinction) from DR8 for secondary
SEC_DR8_I	Dered I magnitude (model - extinction) from DR8 for secondary
SEC_DR8_Z	Dered Z magnitude (model - extinction) from DR8 for secondary
SEC_DR8_SPECZ	Spectral redshift from DR8 for secondary
SEC_DR8_PHOTOZ	Photometric redshift from DR8 for secondary
SEC_DR8_PHOTOZ2	Alternate photometric redshift from DR8 for secondary
