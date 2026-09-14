"""Recompute the submitted benchmark using portable, packaged inputs."""
from pathlib import Path
import argparse,csv,hashlib,json,sys,math

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[1]/'pipeline'))
import benchmark_ranking_metrics as b

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    if args.output_dir.exists():
        raise SystemExit('Use a new output directory to preserve existing results.')
    config=json.loads((ROOT/'evaluation_configuration.json').read_text())
    for item in json.loads((ROOT/'files.sha256.json').read_text()):
        if hashlib.sha256((ROOT/item['file']).read_bytes()).hexdigest()!=item['sha256']:
            raise SystemExit('Evidence hash mismatch: '+item['file'])
    if hashlib.sha256(Path(b.__file__).read_bytes()).hexdigest()!=config['benchmark_code_sha256']:
        raise SystemExit('Benchmark code hash differs from the recorded version.')
    truth=b.load_ground_truth(ROOT/'labels_181.xlsx')
    assert len(truth.allergens)==181 and truth.positive_pair_count==133
    sets={}; all_summary=[]
    for tool,label in [('cross-react','crossreact'),('surface-id','surfaceid')]:
        rankings=b.collect_rankings(sorted((ROOT/(label+'_rankings')).glob('*.csv')),tool)
        assert set(rankings)==set(config['queries']) and len(rankings)==30
        per=[]
        for query,r in sorted(rankings.items()):
            candidates=set(r.targets)
            assert len(r.targets)==len(candidates)==180
            assert query not in candidates and candidates|{query}==set(truth.allergens)
            assert not set(truth.relevant[query])-candidates
            assert sets.setdefault(query,candidates)==candidates
            per.extend(b.compute_query_rows(tool,b.score_name_for_tool(tool),query,r.ranks,r.targets,r.scores,truth.relevant[query],config['ks']))
        summary=b.aggregate_rows(per,config['ks'],tool,b.score_name_for_tool(tool))
        for name,actual in [('per_query_metrics.csv',per),('summary_metrics.csv',summary)]:
            with (ROOT/('expected_'+label)/name).open(encoding='utf-8-sig',newline='') as f:
                expected=list(csv.DictReader(f))
            assert len(actual)==len(expected)
            for left,right in zip(actual,expected):
                for key,value in right.items():
                    try:
                        assert math.isclose(float(left[key]),float(value),rel_tol=0,abs_tol=1e-9),(label,name,key,left[key],value)
                    except ValueError:
                        assert str(left[key])==value,(label,name,key)
        output=args.output_dir/label;output.mkdir(parents=True)
        b.write_csv_atomic(output/'per_query_metrics.csv',b.PER_QUERY_COLUMNS,per)
        b.write_csv_atomic(output/'summary_metrics.csv',b.SUMMARY_COLUMNS,summary)
        all_summary.extend(summary)
    b.write_csv_atomic(args.output_dir/'comparison_summary.csv',b.SUMMARY_COLUMNS,all_summary)
    result={'queries':30,'candidates_per_query':180,'allergens':181,'positive_pairs':133,'original_metrics_match':True,'input_hashes_match':True,'missing_positive_targets':0,'scope':'Saved rankings recomputed; no structural search rerun.'}
    (args.output_dir/'verification.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result))

if __name__=='__main__':main()
