import pandas as pd
import operator as op

def plotLogo(df: pd.DataFrame, output, base_url, top_target_genes: int = 3):
    """
    :param df:
    :param base_url:
    """
    COLUMN_NAME_LOGO = "MotifLogo"
    COLUMN_NAME_MOTIF_ID = "MotifID"
    COLUMN_NAME_TARGETS = "TargetGenes"
    # Make sure the original dataframe is not altered.
    df = df.copy()
    # Add column with URLs to sequence logo.
    def create_url(motif_id):
        return '<img src="{}{}.png" style="max-height:124px;"></img>'.format(base_url, motif_id)
    df[("Enrichment", COLUMN_NAME_LOGO)] = list(map(create_url, df.index.get_level_values(COLUMN_NAME_MOTIF_ID)))
    # Truncate TargetGenes.
    def truncate(col_val):
        return sorted(col_val, key=op.itemgetter(1))[:top_target_genes]
    df[("Enrichment", COLUMN_NAME_TARGETS)] = list(map(truncate, df[("Enrichment", COLUMN_NAME_TARGETS)]))

    html_content = df.to_html(escape=False)
    with open(output, 'w', encoding='utf-8') as html_file:
        html_file.write(html_content)
