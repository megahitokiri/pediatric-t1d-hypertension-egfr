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
    l.Patient_ID,
    l.ITEM_CODE,
    l.ITEM,
    l.OBS_VALUE,
    l.RESULT_DTM,
    TRY_CONVERT(datetime, l.RESULT_DTM) AS ResultDateTime
FROM dbo.Labs_IHC_UUHSC l
JOIN cohort c
  ON c.Patient_ID = l.Patient_ID
WHERE l.OBS_VALUE IS NOT NULL
  AND l.ITEM_CODE IN ('2160-0', '20025')
ORDER BY l.Patient_ID, TRY_CONVERT(datetime, l.RESULT_DTM);
