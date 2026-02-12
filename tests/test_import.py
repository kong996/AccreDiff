def test_import():
    import accrediff as ad
    assert hasattr(ad, "__version__")