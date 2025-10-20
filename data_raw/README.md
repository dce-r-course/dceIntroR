This folder contains two datasets that are used for this course:

1.  `basic_data.rds`
2.  `diag_data.rds`

`basic_data.rds` contains information on

-   `id` (unique personal identifier)

-   `sex` {1 = male, 2 = female}

-   `death_date` (date of death)

-   `birth_date` (date of birth)

`diag_data.rds` contains information on

-   `id` (unique personal identifier)

-   `contact_type` {1 = inpatient, 2 = outpatient}

-   `adm_date` (date of admission)

-   `dis_date` (date of discharge)

-   `diag_a` (primary diagnosis)

-   `diag_b` (secondary diagnosis)

These datasets are thus simulated data aimed at mimicking the Danish registries.

`diag_data.rds` is designed for the disease myocardial infarction of interest. Thus, the `diag_a` and `diag_b` are constructed to reflect a realistic reflection on comorbidity prevalence among myocardial infarction patients. `death_date` is likewise constructed according to a mortality target among patients with myocardial infarction. The dataset is thus not suitable for any other study cohorts.

In `diag_data.rds`, transfers between departments within the same hospital contacts have been taken into account (which will not be the case in *e.g.,* the raw data from the Danish National Patient Registry).

`diag_data.rds` contains information on the following diseases:

| Comorbidity | ICD-10 code |
|----|----|
| Myocardial infarction | DI21 |
| Hypertension | DI10, DI11, DI12, DI13, DI14, DI15 |
| Hypercholesterolaemia | DE78 |
| Diabetes | DE10, DE11, DE12, DE13, DE14 |
| Heart failure | DI50 |
| Random | DX01, DX02, DX03, DX04, DX05, DX06, DX07, DX08, DX09 |
