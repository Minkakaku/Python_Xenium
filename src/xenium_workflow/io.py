from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd


def read_sample_info(path):
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _read_xenium_record(record, data_dir, xenium_jobs, table_key):
    from spatialdata_io import xenium

    raw_dir = Path(data_dir) / record.raw
    sdata = xenium(raw_dir, cells_as_circles=True, n_jobs=xenium_jobs)
    table = sdata.tables[table_key].copy()
    if "cell_id" not in table.obs:
        table.obs["cell_id"] = table.obs_names.astype(str)
    table.obs["cell_id"] = table.obs["cell_id"].astype(str)
    table.obs["fov"] = record.fov
    table.obs["sample"] = record.sample
    table.obs["group"] = record.group
    table.obs_names = pd.Index(
        record.fov + ":" + table.obs["cell_id"].astype(str),
        name="cell_key",
    )
    return table


def _read_roi_map(sample_info, roi_dir):
    rows = []
    for record in sample_info.itertuples(index=False):
        if ".csv_" not in record.sample:
            continue
        tokens = record.sample.split(",")
        groups = record.group.split(",")
        if len(tokens) != len(groups):
            raise ValueError(f"ROI sample/group count differs for FOV {record.fov}")
        for token, group in zip(tokens, groups):
            roi_name, sample = token.split(".csv_", 1)
            roi = pd.read_csv(Path(roi_dir) / f"{roi_name}.csv", skiprows=2)
            roi = roi.rename(columns={"Cell ID": "cell_id"})
            roi["cell_id"] = roi["cell_id"].astype(str)
            roi["fov"] = record.fov
            roi["sample"] = sample
            roi["group"] = group
            rows.append(roi[["cell_id", "fov", "sample", "group"]])
    return pd.concat(rows, ignore_index=True) if rows else None


def _apply_roi_assignments(adata, sample_info, roi_dir):
    roi_map = _read_roi_map(sample_info, roi_dir)
    if roi_map is None:
        return adata

    multi_fovs = sample_info.loc[
        sample_info["sample"].str.contains(".csv_", regex=False), "fov"
    ].tolist()
    roi_map["cell_key"] = roi_map["fov"] + ":" + roi_map["cell_id"]
    roi_map = roi_map.drop_duplicates("cell_key").set_index("cell_key")

    keep = ~adata.obs["fov"].isin(multi_fovs) | adata.obs_names.isin(roi_map.index)
    adata = adata[keep].copy()
    matched = adata.obs_names.intersection(roi_map.index)
    adata.obs.loc[matched, "sample"] = roi_map.loc[matched, "sample"].to_numpy()
    adata.obs.loc[matched, "group"] = roi_map.loc[matched, "group"].to_numpy()
    return adata


def load_xenium_data(
    data_dir,
    sample_info,
    roi_dir=None,
    apply_roi=False,
    workers=1,
    xenium_jobs=6,
    table_key="table",
):
    import anndata

    info = read_sample_info(sample_info)
    records = list(info.itertuples(index=False))

    def read_one(record):
        return _read_xenium_record(record, data_dir, xenium_jobs, table_key)

    if workers > 1:
        with ThreadPoolExecutor(max_workers=workers) as executor:
            tables = list(executor.map(read_one, records))
    else:
        tables = [read_one(record) for record in records]

    adata = anndata.concat(tables, join="inner", merge="same", index_unique=None)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    if apply_roi:
        adata = _apply_roi_assignments(adata, info, roi_dir)
    adata.obs["cell_id_sample"] = (
        adata.obs["cell_id"].astype(str) + "-" + adata.obs["sample"].astype(str)
    )
    return adata
