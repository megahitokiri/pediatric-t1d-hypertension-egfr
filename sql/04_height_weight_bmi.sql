{{COHORT_CTES}},
cohort AS (
    SELECT *
    FROM FINAL
    WHERE AgeAtFirstT1 IS NOT NULL
      AND AgeAtFirstT1 < {{MAX_FIRST_T1_AGE}}
      AND AgeAtLastVisit IS NOT NULL
      AND AgeAtLastVisit < {{MAX_LAST_VISIT_AGE}}
      AND YearsFromFirstT1ToLast IS NOT NULL
      AND YearsFromFirstT1ToLast >= {{MIN_FOLLOWUP_YEARS}}
      AND VisitCount >= {{MIN_VISITS}}
      AND FinalType IN ('T1D', 'IMPUTED_T1D')
)
SELECT
    h.Patient_ID,
    h.UPDB_DIST_ID,
    h.OBSERVATION_DATE,
    h.HEIGHT_CM,
    h.WEIGHT_KG,
    h.BMI
FROM dbo.HEIGHT_WEIGHT_BMI_IHC_UUHSC h
JOIN cohort c
  ON c.Patient_ID = h.Patient_ID
WHERE h.OBSERVATION_DATE IS NOT NULL
ORDER BY h.Patient_ID, h.OBSERVATION_DATE;
