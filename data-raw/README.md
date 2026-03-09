## data-raw content

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

## Simulated Data Design

These datasets are **simulated** to resemble data from Danish health registries. They are specifically tailored to represent patients with **myocardial infarction (MI)**:

-   Diagnoses (`diag_a`, `diag_b`) reflect realistic comorbidity patterns among MI patients.

-   Mortality (`death_date`) is simulated to match expected outcomes for this patient group.

-   The data is **not suitable** for analyses involving other disease cohorts.

Additionally, `diag_data.rds` accounts for **intra-hospital transfers**, meaning multiple department contacts within the same hospital stay are consolidated — unlike raw data from the Danish National Patient Registry.

In `diag_data.rds`, transfers between departments within the same hospital contacts have been taken into account (which will not be the case in *e.g.,* the raw data from the Danish National Patient Registry).

## Diagnoses Included

`diag_data.rds` contains information on the following diseases:

| Comorbidity | ICD-10 code |
|----|----|
| Myocardial infarction | DI21 |
| Hypertension | DI10, DI11, DI12, DI13, DI14, DI15 |
| Hypercholesterolaemia | DE78 |
| Diabetes | DE10, DE11, DE12, DE13, DE14 |
| Heart failure | DI50 |
| Random | DX01, DX02, DX03, DX04, DX05, DX06, DX07, DX08, DX09 |
