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
    bp.Patient_ID,
    bp.OBSERVATION_DATE,
    bp.BP_SYS AS Systolic,
    bp.BP_DIAS AS Diastolic
FROM dbo.BP_IHC_UUHSC bp
JOIN cohort c
  ON c.Patient_ID = bp.Patient_ID
WHERE bp.OBSERVATION_DATE IS NOT NULL
  AND bp.BP_SYS IS NOT NULL
  AND bp.BP_DIAS IS NOT NULL
ORDER BY bp.Patient_ID, bp.OBSERVATION_DATE;
