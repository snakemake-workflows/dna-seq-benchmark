#!/usr/bin/env python3
"""
Tests for calc-vaf.py: VAF calculation from AD field and from two numeric fields.

Expected VAFs for AD-based input (Alt / (Ref + Alt)):
  chr1:100  AD=30,15 -> VAF = 15/45 = 0.3333...
  chr1:200  AD=20,10 -> VAF = 10/30 = 0.3333...
  chr1:300  AD=5,0   -> VAF = 0/5  = 0.0
  chr1:400  AD=0,40  -> VAF = 40/40 = 1.0
  chr1:500  AD=50,0  -> VAF = 0/50 = 0.0

Expected VAFs from fields (NUM / NUM_den):
  chr1:100  15/45 = 0.3333...
  chr1:200  10/30 = 0.3333...
  chr1:300  0/5   = 0.0
  chr1:400  40/40 = 1.0
  chr1:500  0/50  = 0.0
"""

import sys
import os
import numpy as np
from cyvcf2 import VCF

# Load calc-vaf module via importlib (hyphen in filename)
scripts_dir = os.path.join(os.path.dirname(__file__), "..", "workflow", "scripts")
sys.path.insert(0, scripts_dir)
import importlib.util
spec = importlib.util.spec_from_file_location("calc_vaf", os.path.join(scripts_dir, "calc-vaf.py"))
calc_vaf = importlib.util.module_from_spec(spec)
spec.loader.exec_module(calc_vaf)


def test_ad_vaf():
    """Test VAF calculation from AD FORMAT field."""
    vcf_path = os.path.join(os.path.dirname(__file__), "test-vaf-ad.vcf")
    with VCF(vcf_path) as vcf:
        samples = vcf.samples
        expected = {
            "100": 15.0 / 45.0,
            "200": 10.0 / 30.0,
            "300": 0.0 / 5.0,
            "400": 40.0 / 40.0,
            "500": 0.0 / 50.0,
        }
        passed = 0
        failed = 0
        for variant in vcf:
            vaf_array = calc_vaf.calculate_vaf_from_ad(variant, samples)
            assert vaf_array is not None, f"VAF returned None for {variant.CHROM}:{variant.POS}"
            assert vaf_array.shape == (1, 1), f"Wrong shape: {vaf_array.shape}"
            vaf = vaf_array[0, 0]
            exp = expected[str(variant.POS)]
            if np.isclose(vaf, exp, rtol=1e-4):
                print(f"  PASS {variant.CHROM}:{variant.POS} VAF={vaf:.4f} (expected {exp:.4f})")
                passed += 1
            else:
                print(f"  FAIL {variant.CHROM}:{variant.POS} VAF={vaf:.4f} (expected {exp:.4f})")
                failed += 1
    print(f"AD VAF test: {passed} passed, {failed} failed")
    assert failed == 0, f"AD VAF test failed with {failed} errors"


def test_fields_vaf():
    """Test VAF calculation from two FORMAT fields (NUM / NUM_den)."""
    vcf_path = os.path.join(os.path.dirname(__file__), "test-vaf-fields.vcf")
    with VCF(vcf_path) as vcf:
        samples = vcf.samples
        expected = {
            "100": 15.0 / 45.0,
            "200": 10.0 / 30.0,
            "300": 0.0 / 5.0,
            "400": 40.0 / 40.0,
            "500": 0.0 / 50.0,
        }
        passed = 0
        failed = 0
        for variant in vcf:
            vaf_array = calc_vaf.calculate_vaf_from_fields(variant, "FORMAT", "NUM", samples)
            assert vaf_array is not None, f"VAF returned None for {variant.CHROM}:{variant.POS}"
            assert vaf_array.shape == (1, 1), f"Wrong shape: {vaf_array.shape}"
            vaf = vaf_array[0, 0]
            exp = expected[str(variant.POS)]
            if np.isclose(vaf, exp, rtol=1e-4):
                print(f"  PASS {variant.CHROM}:{variant.POS} VAF={vaf:.4f} (expected {exp:.4f})")
                passed += 1
            else:
                print(f"  FAIL {variant.CHROM}:{variant.POS} VAF={vaf:.4f} (expected {exp:.4f})")
                failed += 1
    print(f"Fields VAF test: {passed} passed, {failed} failed")
    assert failed == 0, f"Fields VAF test failed with {failed} errors"


def test_missing_ad():
    """Test that missing AD field returns None."""
    vcf_path = os.path.join(os.path.dirname(__file__), "test-vaf-no-ad.vcf")
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n")
        f.write("chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1\n")
    with VCF(vcf_path) as vcf:
        samples = vcf.samples
        for variant in vcf:
            vaf_array = calc_vaf.calculate_vaf_from_ad(variant, samples)
            assert vaf_array is None, f"Expected None for missing AD, got {vaf_array}"
    os.remove(vcf_path)
    print("  PASS Missing AD field handled correctly")


def test_multi_allele():
    """Test VAF calculation with multi-allelic variant."""
    vcf_path = os.path.join(os.path.dirname(__file__), "test-vaf-multi.vcf")
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String>\n")
        f.write("##FORMAT=<ID=AD,Number=R,Type=Integer>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n")
        # AD = [ref, alt1, alt2] = [50, 30, 20], total=100
        f.write("chr1\t100\t.\tA\tG,T\t.\tPASS\t.\tGT:AD\t0/1:50,30,20\n")
    with VCF(vcf_path) as vcf:
        samples = vcf.samples
        for variant in vcf:
            vaf_array = calc_vaf.calculate_vaf_from_ad(variant, samples)
            assert vaf_array is not None
            assert vaf_array.shape == (1, 2), f"Wrong shape: {vaf_array.shape}"
            assert np.isclose(vaf_array[0, 0], 0.3, rtol=1e-4), f"VAF1 = {vaf_array[0,0]} (expected 0.3)"
            assert np.isclose(vaf_array[0, 1], 0.2, rtol=1e-4), f"VAF2 = {vaf_array[0,1]} (expected 0.2)"
            print(f"  PASS Multi-allele: VAF1={vaf_array[0,0]:.4f}, VAF2={vaf_array[0,1]:.4f}")
    os.remove(vcf_path)
    print("Multi-allele test: PASSED")


def test_info_fields_vaf():
    """Test VAF calculation from two INFO fields."""
    vcf_path = os.path.join(os.path.dirname(__file__), "test-vaf-info.vcf")
    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##FORMAT=<ID=GT,Number=1,Type=String>\n")
        f.write("##INFO=<ID=NUM,Number=1,Type=Integer>\n")
        f.write("##INFO=<ID=NUM_den,Number=1,Type=Integer>\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample1\n")
        f.write("chr1\t100\t.\tA\tG\t.\tPASS\tNUM=15;NUM_den=45\tGT\t0/1\n")
    with VCF(vcf_path) as vcf:
        samples = vcf.samples
        for variant in vcf:
            vaf_array = calc_vaf.calculate_vaf_from_fields(variant, "INFO", "NUM", samples)
            assert vaf_array is not None
            assert vaf_array.shape == (1, 1), f"Wrong shape: {vaf_array.shape}"
            assert np.isclose(vaf_array[0, 0], 15.0 / 45.0, rtol=1e-4)
            print(f"  PASS INFO fields: VAF={vaf_array[0,0]:.4f}")
    os.remove(vcf_path)
    print("INFO fields test: PASSED")


def test_add_vaf_to_vcf_ad():
    """Test the full add_vaf_to_vcf function with AD-based calculation."""
    import tempfile

    input_path = os.path.join(os.path.dirname(__file__), "test-vaf-ad.vcf")
    with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as f:
        output_path = f.name

    try:
        calc_vaf.add_vaf_to_vcf(input_path, output_path, calculate_from_ad=True)

        # Verify output has VAF field
        with VCF(output_path) as vcf:
            samples = vcf.samples
            # Check header has VAF format by trying to read it from a variant
            # (cyvcf2 doesn't expose header_lines directly)
            var = vcf.__next__()
            try:
                var.format("VAF")  # Should not raise
            except KeyError:
                assert False, "VAF not in header"

            # Verify VAF values
            for variant in vcf:
                vaf = variant.format("VAF")
                assert vaf is not None, f"VAF missing for {variant.CHROM}:{variant.POS}"
                assert vaf.shape == (1, 1), f"VAF wrong shape: {vaf.shape}"
                print(f"  PASS {variant.CHROM}:{variant.POS} output VAF={vaf[0,0]:.4f}")
        print("Full add_vaf_to_vcf (AD) test: PASSED")
    finally:
        if os.path.exists(output_path):
            os.remove(output_path)


def test_add_vaf_to_vcf_fields():
    """Test the full add_vaf_to_vcf function with custom field calculation."""
    import tempfile

    input_path = os.path.join(os.path.dirname(__file__), "test-vaf-fields.vcf")
    with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as f:
        output_path = f.name

    try:
        calc_vaf.add_vaf_to_vcf(
            input_path, output_path,
            calculate_from_ad=False,
            num_field="FORMAT", num_name="NUM",
            den_field="FORMAT", den_name="NUM_den",
        )

        # Verify output has VAF field
        with VCF(output_path) as vcf:
            var = vcf.__next__()
            try:
                var.format("VAF")  # Should not raise
            except KeyError:
                assert False, "VAF not in header"

            for variant in vcf:
                vaf = variant.format("VAF")
                assert vaf is not None, f"VAF missing for {variant.CHROM}:{variant.POS}"
                assert vaf.shape == (1, 1), f"VAF wrong shape: {vaf.shape}"
                print(f"  PASS {variant.CHROM}:{variant.POS} output VAF={vaf[0,0]:.4f}")
        print("Full add_vaf_to_vcf (fields) test: PASSED")
    finally:
        if os.path.exists(output_path):
            os.remove(output_path)


if __name__ == "__main__":
    os.chdir(os.path.dirname(__file__))

    tests = [
        ("AD VAF", test_ad_vaf),
        ("Fields VAF", test_fields_vaf),
        ("Missing AD", test_missing_ad),
        ("Multi-allele", test_multi_allele),
        ("INFO fields", test_info_fields_vaf),
        ("Full add_vaf_to_vcf (AD)", test_add_vaf_to_vcf_ad),
        ("Full add_vaf_to_vcf (fields)", test_add_vaf_to_vcf_fields),
    ]

    for name, test_func in tests:
        print("=" * 60)
        print(f"Testing {name}")
        print("=" * 60)
        try:
            test_func()
            print(f"✓ {name} PASSED")
        except Exception as e:
            print(f"✗ {name} FAILED: {e}")
            sys.exit(1)
        print()

    print("=" * 60)
    print("All tests PASSED!")
    print("=" * 60)
