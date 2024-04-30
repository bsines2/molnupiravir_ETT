

@transform_pandas(
    Output(rid="ri.foundry.main.dataset.b9336838-b304-462e-911b-331bf5f44ab7"),
    COVID_Patient_Summary_Table=Input(rid="ri.foundry.main.dataset.41ea332f-c59a-4b3d-9ccf-529854028b78")
)
SELECT COUNT(DISTINCT person_id) AS N,
shift_date_yn
FROM COVID_Patient_Summary_Table
GROUP BY 2

