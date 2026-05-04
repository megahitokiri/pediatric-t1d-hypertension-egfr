WITH DX AS (
    SELECT
        p.Patient_ID,
        p.SEX_CD,
        p.BIRTH_DT,
        d.ICD_CODE,
        d.TYPE_1OR2,
        d.ADM_DATE,
        d.AGE_AT_ADM,
        CASE
            WHEN d.TYPE_1OR2 = '1'
                 OR d.ICD_CODE LIKE 'E10%'
                 OR (d.ICD_CODE LIKE '250.%' AND RIGHT(d.ICD_CODE, 1) IN ('1', '3'))
            THEN 1 ELSE 0
        END AS is_t1,
        CASE
            WHEN d.TYPE_1OR2 = '2'
                 OR d.ICD_CODE LIKE 'E11%'
                 OR (d.ICD_CODE LIKE '250.%' AND RIGHT(d.ICD_CODE, 1) IN ('0', '2'))
            THEN 1 ELSE 0
        END AS is_t2
    FROM dbo.Patients_IHC_UUHSC p
    JOIN dbo.Diagnoses_IHC_UUHSC d
      ON p.Patient_ID = d.Patient_ID
    WHERE d.ADM_DATE IS NOT NULL
),
AGG AS (
    SELECT
        Patient_ID,
        MAX(SEX_CD) AS SEX_CD,
        MIN(BIRTH_DT) AS BIRTH_DT,
        SUM(is_t1) AS Type1_count,
        SUM(is_t2) AS Type2_count,
        MIN(CASE WHEN is_t1 = 1 THEN ADM_DATE END) AS first_t1_date,
        MIN(CASE WHEN is_t1 = 1 THEN AGE_AT_ADM END) AS first_t1_age_stored,
        MAX(ADM_DATE) AS last_visit_date,
        COUNT(DISTINCT ADM_DATE) AS VisitCount
    FROM DX
    GROUP BY Patient_ID
),
LAST_AGE AS (
    SELECT
        a.*,
        x.AGE_AT_ADM AS last_age_stored
    FROM AGG a
    OUTER APPLY (
        SELECT TOP 1 dx.AGE_AT_ADM
        FROM DX dx
        WHERE dx.Patient_ID = a.Patient_ID
          AND dx.ADM_DATE = a.last_visit_date
          AND dx.AGE_AT_ADM IS NOT NULL
        ORDER BY dx.AGE_AT_ADM DESC
    ) x
),
AGE AS (
    SELECT
        *,
        CASE
            WHEN first_t1_age_stored IS NOT NULL THEN first_t1_age_stored
            WHEN BIRTH_DT IS NOT NULL AND first_t1_date IS NOT NULL
                 THEN DATEDIFF(day, BIRTH_DT, first_t1_date) / 365.25
            ELSE NULL
        END AS AgeAtFirstT1,
        CASE
            WHEN last_age_stored IS NOT NULL THEN last_age_stored
            WHEN BIRTH_DT IS NOT NULL AND last_visit_date IS NOT NULL
                 THEN DATEDIFF(day, BIRTH_DT, last_visit_date) / 365.25
            ELSE NULL
        END AS AgeAtLastVisit,
        CASE
            WHEN first_t1_date IS NOT NULL AND last_visit_date IS NOT NULL
                 THEN DATEDIFF(day, first_t1_date, last_visit_date) / 365.25
            ELSE NULL
        END AS YearsFromFirstT1ToLast
    FROM LAST_AGE
),
FINAL AS (
    SELECT
        Patient_ID,
        SEX_CD,
        BIRTH_DT,
        Type1_count,
        Type2_count,
        first_t1_date,
        last_visit_date,
        AgeAtFirstT1,
        AgeAtLastVisit,
        VisitCount,
        YearsFromFirstT1ToLast,
        CAST(Type1_count AS float) / NULLIF(CAST(Type2_count AS float), 0) AS Type1_Type2_Ratio,
        CASE
            WHEN Type1_count >= 1 AND (Type2_count IS NULL OR Type2_count = 0) THEN 'T1D'
            WHEN Type2_count >= 1 AND (Type1_count IS NULL OR Type1_count = 0) THEN 'T2D'
            WHEN (CAST(Type1_count AS float) / NULLIF(CAST(Type2_count AS float), 0)) BETWEEN 0.33 AND 3.33 THEN 'T3D'
            WHEN (CAST(Type1_count AS float) / NULLIF(CAST(Type2_count AS float), 0)) > 3.33 THEN 'IMPUTED_T1D'
            WHEN (CAST(Type1_count AS float) / NULLIF(CAST(Type2_count AS float), 0)) < 0.33 THEN 'IMPUTED_T2D'
            ELSE 'OTHER'
        END AS FinalType
    FROM AGE
)
