#!/usr/bin/env python3
"""
PACE Comprehensive Test Suite
测试所有模块，包括机器学习部分

Author: Linyong Shen @ Northwest A&F University
"""

import os
import sys
import tempfile
import shutil
import traceback
from pathlib import Path

# Get PACE root directory (where this script is located)
PACE_ROOT = Path(__file__).parent.absolute()
print(f"PACE_ROOT: {PACE_ROOT}")

# Add paths for module imports
sys.path.insert(0, str(PACE_ROOT / "workflow" / "scripts"))
sys.path.insert(0, str(PACE_ROOT / "scripts"))

# Test results tracking
PASSED = []
FAILED = []

def test(name):
    """Decorator for test functions"""
    def decorator(func):
        def wrapper(*args, **kwargs):
            try:
                result = func(*args, **kwargs)
                PASSED.append(name)
                print(f"  ✓ PASSED: {name}")
                return result
            except Exception as e:
                FAILED.append((name, str(e)))
                print(f"  ✗ FAILED: {name}")
                print(f"    Error: {e}")
                traceback.print_exc()
                return None
        return wrapper
    return decorator


class TestData:
    """Create test data for all modules"""
    
    def __init__(self, test_dir):
        self.test_dir = Path(test_dir)
        self.test_dir.mkdir(exist_ok=True)
        
    def create_all(self):
        """Create all test data files"""
        self.create_chrom_sizes()
        self.create_peaks()
        self.create_genes()
        self.create_tss()
        self.create_blacklist()
        self.create_expression()
        self.create_enhancers()
        self.create_gene_list()
        self.create_accessibility()
        self.create_predictions()
        self.create_validation()
        return self
        
    def create_chrom_sizes(self):
        with open(self.test_dir / "chrom.sizes", 'w') as f:
            f.write("chr1\t100000\nchr2\t80000\nchr3\t60000\n")
    
    def create_peaks(self):
        with open(self.test_dir / "test_peaks.narrowPeak", 'w') as f:
            # chr start end name score strand signalValue pValue qValue peak(summit offset)
            f.write("chr1\t1000\t2000\tpeak1\t500\t.\t10.5\t5.0\t3.0\t500\n")
            f.write("chr1\t5000\t6000\tpeak2\t800\t.\t15.2\t8.0\t5.0\t500\n")
            f.write("chr2\t2000\t3000\tpeak3\t600\t.\t12.0\t6.0\t4.0\t500\n")
            f.write("chr2\t10000\t11000\tpeak4\t400\t.\t8.5\t4.0\t2.5\t500\n")
            f.write("chr3\t3000\t4000\tpeak5\t700\t.\t14.0\t7.0\t4.5\t500\n")
    
    def create_genes(self):
        with open(self.test_dir / "test_genes.bed", 'w') as f:
            f.write("chr1\t8000\t15000\tGene1\t0\t+\n")
            f.write("chr2\t5000\t12000\tGene2\t0\t-\n")
            f.write("chr3\t1000\t8000\tGene3\t0\t+\n")
    
    def create_tss(self):
        with open(self.test_dir / "test_tss.bed", 'w') as f:
            f.write("chr1\t7500\t8500\tGene1_TSS\t0\t+\n")
            f.write("chr2\t11500\t12500\tGene2_TSS\t0\t-\n")
            f.write("chr3\t500\t1500\tGene3_TSS\t0\t+\n")
    
    def create_blacklist(self):
        with open(self.test_dir / "test_blacklist.bed", 'w') as f:
            f.write("chr1\t50000\t51000\tblacklist1\n")
    
    def create_expression(self):
        with open(self.test_dir / "test_expression.tsv", 'w') as f:
            f.write("gene_id\tgene_name\tTPM\n")
            f.write("ENSG001\tGene1\t50.5\n")
            f.write("ENSG002\tGene2\t120.3\n")
            f.write("ENSG003\tGene3\t0.5\n")
    
    def create_enhancers(self):
        with open(self.test_dir / "test_enhancers.txt", 'w') as f:
            f.write("chr\tstart\tend\tname\tactivity_base\n")
            f.write("chr1\t1000\t1500\tenh1\t1.5\n")
            f.write("chr1\t5000\t5500\tenh2\t2.0\n")
            f.write("chr2\t2000\t2500\tenh3\t1.8\n")
            f.write("chr2\t10000\t10500\tenh4\t1.2\n")
            f.write("chr3\t3000\t3500\tenh5\t1.6\n")
    
    def create_gene_list(self):
        with open(self.test_dir / "test_gene_list.txt", 'w') as f:
            f.write("chr\tstart\tend\tname\tExpression\tTSS\tstrand\n")
            f.write("chr1\t8000\t15000\tGene1\t50.5\t8000\t+\n")
            f.write("chr2\t5000\t12000\tGene2\t120.3\t12000\t-\n")
            f.write("chr3\t1000\t8000\tGene3\t0.5\t1000\t+\n")
    
    def create_accessibility(self):
        with open(self.test_dir / "test_accessibility.tagAlign", 'w') as f:
            # Create reads around peaks
            for pos in [1200, 1300, 1400, 5200, 5300, 5400]:
                f.write(f"chr1\t{pos}\t{pos+50}\tN\t1000\t+\n")
            for pos in [2200, 2300, 10200, 10300]:
                f.write(f"chr2\t{pos}\t{pos+50}\tN\t1000\t+\n")
            for pos in [3200, 3300, 3400, 3500]:
                f.write(f"chr3\t{pos}\t{pos+50}\tN\t1000\t+\n")
    
    def create_predictions(self):
        """Create predictions file for ML testing"""
        with open(self.test_dir / "test_predictions.tsv", 'w') as f:
            f.write("chr\tstart\tend\tname\tTargetGene\tTargetGeneTSS\tdistance\t")
            f.write("activity\tcontact\tABC.Score\tH3K27ac\tATAC\n")
            # Generate 100 predictions
            import random
            random.seed(42)
            for i in range(100):
                chrom = f"chr{(i % 3) + 1}"
                start = random.randint(1000, 50000)
                end = start + 500
                name = f"E{i}"
                gene = f"Gene{(i % 10) + 1}"
                tss = random.randint(1000, 50000)
                distance = abs(start - tss)
                activity = random.uniform(0.1, 2.0)
                contact = random.uniform(0.1, 0.8)
                abc_score = random.uniform(0.01, 0.3)
                h3k27ac = random.uniform(0, 1)
                atac = random.uniform(0, 1)
                f.write(f"{chrom}\t{start}\t{end}\t{name}\t{gene}\t{tss}\t{distance}\t")
                f.write(f"{activity:.4f}\t{contact:.4f}\t{abc_score:.4f}\t{h3k27ac:.4f}\t{atac:.4f}\n")
    
    def create_validation(self):
        """Create validation data for ML testing"""
        with open(self.test_dir / "test_validation.tsv", 'w') as f:
            f.write("enhancer\tgene\tvalidated\n")
            # Read predictions and mark some as validated
            import random
            random.seed(42)
            with open(self.test_dir / "test_predictions.tsv") as pred_f:
                next(pred_f)  # Skip header
                for i, line in enumerate(pred_f):
                    if i >= 50:
                        break
                    parts = line.strip().split('\t')
                    chrom, start, end = parts[0], parts[1], parts[2]
                    gene = parts[4]
                    enhancer = f"{chrom}:{start}-{end}"
                    # Higher ABC score = more likely validated
                    abc_score = float(parts[9])
                    validated = 1 if abc_score > 0.15 else 0
                    f.write(f"{enhancer}\t{gene}\t{validated}\n")


def run_module_tests():
    """Test core modules"""
    print("\n" + "=" * 60)
    print("测试核心模块导入")
    print("=" * 60)
    
    @test("import tools")
    def test_import_tools():
        import tools
        return tools
    
    @test("import peaks")
    def test_import_peaks():
        import peaks
        return peaks
    
    @test("import neighborhoods")
    def test_import_neighborhoods():
        import neighborhoods
        return neighborhoods
    
    @test("import hic")
    def test_import_hic():
        import hic
        return hic
    
    @test("import predictor")
    def test_import_predictor():
        import predictor
        return predictor
    
    @test("import metrics")
    def test_import_metrics():
        import metrics
        return metrics
    
    @test("import ml_integration")
    def test_import_ml():
        import ml_integration
        return ml_integration
    
    test_import_tools()
    test_import_peaks()
    test_import_neighborhoods()
    test_import_hic()
    test_import_predictor()
    test_import_metrics()
    test_import_ml()


def run_function_tests(test_dir):
    """Test core functions"""
    print("\n" + "=" * 60)
    print("测试核心函数")
    print("=" * 60)
    
    import numpy as np
    import pandas as pd
    
    @test("tools.read_bed")
    def test_read_bed():
        from tools import read_bed
        df = read_bed(test_dir / "test_peaks.narrowPeak")
        assert len(df) == 5, f"Expected 5 peaks, got {len(df)}"
        return df
    
    @test("tools.load_chromosome_sizes")
    def test_load_chrom_sizes():
        from tools import load_chromosome_sizes
        sizes = load_chromosome_sizes(test_dir / "chrom.sizes")
        assert 'chr1' in sizes, "chr1 not found"
        assert sizes['chr1'] == 100000, f"Wrong chr1 size: {sizes['chr1']}"
        return sizes
    
    @test("tools.calculate_distance")
    def test_calculate_distance():
        from tools import calculate_distance
        d = calculate_distance(1000, 5000)
        assert d == 4000, f"Expected 4000, got {d}"
        return d
    
    @test("tools.power_law_contact")
    def test_power_law():
        from tools import power_law_contact
        c = power_law_contact(10000)
        assert 0 < c < 1, f"Contact should be 0-1, got {c}"
        return c
    
    @test("tools.geometric_mean")
    def test_geometric_mean():
        from tools import geometric_mean
        gm = geometric_mean([4, 9])
        assert abs(gm - 6.0) < 0.01, f"Expected 6.0, got {gm}"
        return gm
    
    @test("peaks.load_narrowpeak")
    def test_load_narrowpeak():
        from peaks import load_narrowpeak
        df = load_narrowpeak(test_dir / "test_peaks.narrowPeak")
        assert len(df) == 5, f"Expected 5 peaks, got {len(df)}"
        assert 'signalValue' in df.columns, "signalValue column missing"
        return df
    
    @test("peaks.select_strongest_peaks")
    def test_select_peaks():
        from peaks import load_narrowpeak, select_strongest_peaks
        df = load_narrowpeak(test_dir / "test_peaks.narrowPeak")
        selected = select_strongest_peaks(df, n_peaks=3)
        assert len(selected) == 3, f"Expected 3 peaks, got {len(selected)}"
        return selected
    
    @test("peaks.extend_peaks_from_summit")
    def test_extend_peaks():
        from peaks import load_narrowpeak, extend_peaks_from_summit
        df = load_narrowpeak(test_dir / "test_peaks.narrowPeak")
        extended = extend_peaks_from_summit(df, extend_bp=250)
        # Check width is ~500bp
        widths = extended['end'] - extended['start']
        assert all(widths == 500), "Extended peaks should be 500bp wide"
        return extended
    
    @test("hic.ContactEstimator")
    def test_contact_estimator():
        from hic import ContactEstimator
        estimator = ContactEstimator()
        # API: estimate(chrom, pos1, pos2)
        contact = estimator.estimate('chr1', 10000, 20000)
        assert 0 < contact < 1, f"Contact should be 0-1, got {contact}"
        return contact
    
    test_read_bed()
    test_load_chrom_sizes()
    test_calculate_distance()
    test_power_law()
    test_geometric_mean()
    test_load_narrowpeak()
    test_select_peaks()
    test_extend_peaks()
    test_contact_estimator()


def run_cli_tests():
    """Test command-line scripts"""
    print("\n" + "=" * 60)
    print("测试命令行脚本 --help")
    print("=" * 60)
    
    import subprocess
    scripts_dir = PACE_ROOT / "workflow" / "scripts"
    
    scripts = [
        "pace_candidate_regions.py",
        "pace_neighborhoods.py", 
        "pace_predict.py",
        "pace_filter.py",
        "pace_metrics.py"
    ]
    
    for script in scripts:
        @test(f"{script} --help")
        def test_script_help(s=script):
            result = subprocess.run(
                [sys.executable, str(scripts_dir / s), "--help"],
                capture_output=True, text=True
            )
            assert result.returncode == 0, f"Failed: {result.stderr}"
            return result.stdout
        test_script_help()
    
    # Test ML script
    @test("pace_ml.py --help")
    def test_ml_help():
        result = subprocess.run(
            [sys.executable, str(PACE_ROOT / "scripts" / "pace_ml.py"), "--help"],
            capture_output=True, text=True
        )
        assert result.returncode == 0, f"Failed: {result.stderr}"
        return result.stdout
    test_ml_help()


def run_pipeline_tests(test_dir):
    """Test complete pipeline"""
    print("\n" + "=" * 60)
    print("测试完整流程")
    print("=" * 60)
    
    import subprocess
    scripts_dir = PACE_ROOT / "workflow" / "scripts"
    output_dir = test_dir / "output"
    output_dir.mkdir(exist_ok=True)
    
    # Step 1: Candidate regions
    @test("Step 1: pace_candidate_regions.py")
    def test_step1():
        cmd = [
            sys.executable, str(scripts_dir / "pace_candidate_regions.py"),
            "--narrowPeak", str(test_dir / "test_peaks.narrowPeak"),
            "--chrom_sizes", str(test_dir / "chrom.sizes"),
            "--output", str(output_dir / "candidates.bed"),
            "--nStrongestPeaks", "100",
            "--peakExtendFromSummit", "250"
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Failed: {result.stderr}"
        assert (output_dir / "candidates.bed").exists(), "Output not created"
        return result
    
    # Step 2: Neighborhoods
    @test("Step 2: pace_neighborhoods.py")
    def test_step2():
        cmd = [
            sys.executable, str(scripts_dir / "pace_neighborhoods.py"),
            "--candidate_regions", str(output_dir / "candidates.bed"),
            "--genes", str(test_dir / "test_genes.bed"),
            "--chrom_sizes", str(test_dir / "chrom.sizes"),
            "--output_dir", str(output_dir / "neighborhoods"),
            "--accessibility_file", str(test_dir / "test_accessibility.tagAlign")
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Failed: {result.stderr}"
        assert (output_dir / "neighborhoods" / "EnhancerList.txt").exists()
        return result
    
    # Step 3: Predictions
    @test("Step 3: pace_predict.py")
    def test_step3():
        cmd = [
            sys.executable, str(scripts_dir / "pace_predict.py"),
            "--enhancers", str(output_dir / "neighborhoods" / "EnhancerList.txt"),
            "--genes", str(output_dir / "neighborhoods" / "GeneList.txt"),
            "--output", str(output_dir / "predictions.tsv"),
            "--max_distance", "50000",
            "--threshold", "0.01"
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Failed: {result.stderr}"
        assert (output_dir / "predictions.tsv").exists()
        return result
    
    # Step 4: Filter
    @test("Step 4: pace_filter.py")
    def test_step4():
        cmd = [
            sys.executable, str(scripts_dir / "pace_filter.py"),
            "--predictions", str(output_dir / "predictions.tsv"),
            "--output", str(output_dir / "filtered.tsv"),
            "--threshold", "0.02"
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Failed: {result.stderr}"
        return result
    
    # Step 5: Metrics
    @test("Step 5: pace_metrics.py")
    def test_step5():
        cmd = [
            sys.executable, str(scripts_dir / "pace_metrics.py"),
            "--predictions", str(output_dir / "predictions.tsv"),
            "--output_dir", str(output_dir / "qc"),
            "--sample_name", "test"
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Failed: {result.stderr}"
        return result
    
    test_step1()
    test_step2()
    test_step3()
    test_step4()
    test_step5()


def run_ml_tests(test_dir):
    """Test machine learning module"""
    print("\n" + "=" * 60)
    print("测试机器学习模块")
    print("=" * 60)
    
    import subprocess
    ml_script = PACE_ROOT / "scripts" / "pace_ml.py"
    output_dir = test_dir / "ml_output"
    output_dir.mkdir(exist_ok=True)
    
    # Check sklearn availability
    @test("ML: sklearn availability")
    def test_sklearn():
        try:
            import sklearn
            print(f"    sklearn version: {sklearn.__version__}")
            return sklearn.__version__
        except ImportError:
            raise ImportError("scikit-learn not installed")
    
    sklearn_version = test_sklearn()
    
    if sklearn_version is None:
        print("  ⚠ Skipping ML tests - sklearn not available")
        return
    
    # Test ML model classes
    @test("ML: MLConfig dataclass")
    def test_mlconfig():
        from ml_integration import MLConfig
        config = MLConfig(
            model_type="gradient_boosting",
            n_estimators=50,
            max_depth=3
        )
        assert config.model_type == "gradient_boosting"
        return config
    
    @test("ML: PACEMLPredictor init")
    def test_predictor_init():
        from ml_integration import PACEMLPredictor, MLConfig
        config = MLConfig(n_estimators=10, max_depth=2)
        predictor = PACEMLPredictor(config)
        assert predictor.model is None  # Not trained yet
        return predictor
    
    @test("ML: FeatureEngineering")
    def test_feature_eng():
        import pandas as pd
        import numpy as np
        from ml_integration import FeatureEngineering
        
        df = pd.DataFrame({
            'distance': [1000, 5000, 10000],
            'activity': [1.0, 0.5, 0.8],
            'contact': [0.5, 0.3, 0.2],
            'H3K27ac': [0.8, 0.6, 0.4],
            'ATAC': [0.7, 0.5, 0.3]
        })
        
        result = FeatureEngineering.create_features(df)
        assert 'log_distance' in result.columns
        assert 'is_proximal' in result.columns
        return result
    
    @test("ML: Model training (synthetic)")
    def test_model_training():
        import numpy as np
        import pandas as pd
        from ml_integration import PACEMLPredictor, MLConfig
        
        # Create synthetic data
        np.random.seed(42)
        n = 100
        X = np.random.rand(n, 5)
        y = (X[:, 0] + X[:, 1] > 1).astype(int)  # Simple rule
        
        config = MLConfig(n_estimators=10, max_depth=2, n_cv_folds=3)
        predictor = PACEMLPredictor(config)
        
        feature_names = ['f1', 'f2', 'f3', 'f4', 'f5']
        metrics = predictor.train(X, y, feature_names)
        
        assert 'cv_auc_mean' in metrics
        assert metrics['cv_auc_mean'] > 0.5  # Better than random
        print(f"    CV AUC: {metrics['cv_auc_mean']:.3f}")
        return metrics
    
    @test("ML: Model prediction")
    def test_model_prediction():
        import numpy as np
        from ml_integration import PACEMLPredictor, MLConfig
        
        np.random.seed(42)
        n = 100
        X = np.random.rand(n, 5)
        y = (X[:, 0] + X[:, 1] > 1).astype(int)
        
        config = MLConfig(n_estimators=10, max_depth=2, n_cv_folds=3)
        predictor = PACEMLPredictor(config)
        predictor.train(X, y, ['f1', 'f2', 'f3', 'f4', 'f5'])
        
        # Predict
        probs, preds = predictor.predict(X[:10])
        assert len(probs) == 10
        assert all(0 <= p <= 1 for p in probs)
        return probs
    
    @test("ML: Feature importance")
    def test_feature_importance():
        import numpy as np
        from ml_integration import PACEMLPredictor, MLConfig
        
        np.random.seed(42)
        n = 100
        X = np.random.rand(n, 5)
        y = (X[:, 0] + X[:, 1] > 1).astype(int)
        
        config = MLConfig(n_estimators=10, max_depth=2, n_cv_folds=3)
        predictor = PACEMLPredictor(config)
        predictor.train(X, y, ['f1', 'f2', 'f3', 'f4', 'f5'])
        
        importance = predictor.get_feature_importance()
        assert len(importance) == 5
        assert importance['importance'].sum() > 0
        print(f"    Top feature: {importance.iloc[0]['feature']}")
        return importance
    
    test_mlconfig()
    test_predictor_init()
    test_feature_eng()
    test_model_training()
    test_model_prediction()
    test_feature_importance()
    
    # Test CLI ML commands
    @test("ML CLI: train command")
    def test_ml_train():
        cmd = [
            sys.executable, str(ml_script), "train",
            "--predictions", str(test_dir / "test_predictions.tsv"),
            "--validation", str(test_dir / "test_validation.tsv"),
            "--output", str(output_dir / "model.pkl"),
            "--n_estimators", "10",
            "--max_depth", "2",
            "--cv_folds", "3"
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"    STDOUT: {result.stdout}")
            print(f"    STDERR: {result.stderr}")
        assert result.returncode == 0, f"Failed: {result.stderr}"
        assert (output_dir / "model.pkl").exists(), "Model not saved"
        return result
    
    @test("ML CLI: predict command")
    def test_ml_predict():
        cmd = [
            sys.executable, str(ml_script), "predict",
            "--predictions", str(test_dir / "test_predictions.tsv"),
            "--model", str(output_dir / "model.pkl"),
            "--output", str(output_dir / "predictions_ml.tsv"),
            "--abc_weight", "0.5"
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"    STDOUT: {result.stdout}")
            print(f"    STDERR: {result.stderr}")
        assert result.returncode == 0, f"Failed: {result.stderr}"
        assert (output_dir / "predictions_ml.tsv").exists()
        
        # Check output has ML columns
        import pandas as pd
        df = pd.read_csv(output_dir / "predictions_ml.tsv", sep='\t')
        assert 'ML.Score' in df.columns, "ML.Score column missing"
        assert 'Combined.Score' in df.columns, "Combined.Score column missing"
        print(f"    ML.Score range: {df['ML.Score'].min():.4f} - {df['ML.Score'].max():.4f}")
        return df
    
    @test("ML: Random Forest model")
    def test_random_forest():
        cmd = [
            sys.executable, str(ml_script), "train",
            "--predictions", str(test_dir / "test_predictions.tsv"),
            "--validation", str(test_dir / "test_validation.tsv"),
            "--output", str(output_dir / "model_rf.pkl"),
            "--model_type", "random_forest",
            "--n_estimators", "10",
            "--max_depth", "2",
            "--cv_folds", "3"
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        assert result.returncode == 0, f"Failed: {result.stderr}"
        return result
    
    test_ml_train()
    test_ml_predict()
    test_random_forest()


def check_dependencies():
    """Check and report dependency versions"""
    print("\n" + "=" * 60)
    print("检查依赖版本")
    print("=" * 60)
    
    dependencies = [
        ('python', lambda: sys.version.split()[0]),
        ('numpy', lambda: __import__('numpy').__version__),
        ('pandas', lambda: __import__('pandas').__version__),
        ('scipy', lambda: __import__('scipy').__version__),
        ('matplotlib', lambda: __import__('matplotlib').__version__),
        ('seaborn', lambda: __import__('seaborn').__version__),
        ('sklearn', lambda: __import__('sklearn').__version__),
        ('pybedtools', lambda: __import__('pybedtools').__version__),
    ]
    
    versions = {}
    for name, get_version in dependencies:
        try:
            ver = get_version()
            print(f"  ✓ {name}: {ver}")
            versions[name] = ver
        except Exception as e:
            print(f"  ✗ {name}: NOT INSTALLED ({e})")
            versions[name] = None
    
    return versions


def main():
    print("\n" + "=" * 60)
    print("PACE 全面测试套件")
    print("=" * 60)
    
    # Check dependencies
    versions = check_dependencies()
    
    # Create test directory
    test_dir = Path(tempfile.mkdtemp(prefix="pace_test_"))
    print(f"\n测试目录: {test_dir}")
    
    try:
        # Create test data
        print("\n创建测试数据...")
        data = TestData(test_dir).create_all()
        print("  ✓ 测试数据创建完成")
        
        # Run tests
        run_module_tests()
        run_function_tests(test_dir)
        run_cli_tests()
        run_pipeline_tests(test_dir)
        run_ml_tests(test_dir)
        
    finally:
        # Cleanup
        print(f"\n清理测试目录: {test_dir}")
        shutil.rmtree(test_dir, ignore_errors=True)
    
    # Summary
    print("\n" + "=" * 60)
    print("测试总结")
    print("=" * 60)
    print(f"\n总计: {len(PASSED) + len(FAILED)} 测试")
    print(f"通过: {len(PASSED)} ✓")
    print(f"失败: {len(FAILED)} ✗")
    
    if FAILED:
        print("\n失败的测试:")
        for name, error in FAILED:
            print(f"  ✗ {name}")
            print(f"    Error: {error}")
    
    print(f"\n通过率: {100 * len(PASSED) / (len(PASSED) + len(FAILED)):.1f}%")
    
    return len(FAILED) == 0


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
