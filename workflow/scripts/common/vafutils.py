"""Shared VAF extraction utilities."""
import numpy as np


def get_vaf_from_record(record, field, name):
    """Extract VAF value from a record, handling FORMAT and INFO fields.

    Returns a scalar float VAF (first value if array-like), or nan if not found.

    Args:
        record: A pysam VCF record object with .info and .samples.
        field: Either "INFO" or a FORMAT field identifier.
        name: The field ID to look up.

    Returns:
        float: Scalar VAF value, or float('nan') if unavailable.
    """
    try:
        if field == "INFO":
            vaf = record.info.get(name)
        else:
            sample_name = list(record.samples.keys())[0]
            vaf = record.samples[sample_name].get(name)
    except (KeyError, IndexError, AttributeError):
        return float('nan')

    # Handle array-like values (e.g., multi-allelic variants)
    if isinstance(vaf, (list, tuple, np.ndarray, np.generic)):
        if hasattr(vaf, 'item'):
            vaf = vaf.item()
        elif len(vaf) > 0:
            vaf = vaf[0]
        else:
            return float('nan')

    # Convert string percentages (e.g., "50%") to float (e.g., 0.5)
    if isinstance(vaf, str):
        vaf = float(vaf.replace("%", "")) / 100

    # Handle numpy types
    if hasattr(vaf, 'item'):
        vaf = vaf.item()

    try:
        return float(vaf)
    except (ValueError, TypeError):
        return float('nan')
