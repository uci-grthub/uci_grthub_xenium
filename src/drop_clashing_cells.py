def drop_clashing_cells_from_spatial_adata(adata, obs_key='age'):
    """Drop clashing cells from an AnnData object.

    Clashing cells are defined as cells that have the same spatial coordinates.
    This function retains only the first occurrence of each unique coordinate pair.

    Parameters
    ----------
    adata : AnnData
        The input AnnData object containing spatial transcriptomics data.

    Returns
    -------
    AnnData
        A new AnnData object with clashing cells removed.
    """

    # Extract spatial coordinates
    coords = adata.obsm['spatial']

    obs_df = adata.obs[obs_key].to_frame()

    # Create a DataFrame to identify duplicates
    coord_df = pd.DataFrame(coords, columns=['x', 'y'], index = adata.obs.index)

    coord_df = pd.concat([coord_df, obs_df], axis=1)

    # Subset the AnnData object to retain only unique cells
    adata_unique = adata[coord_df.drop_duplicates().index].copy()
    
    return adata_unique