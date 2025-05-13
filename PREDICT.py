import sys
import os
import argparse
sys.path.append('./Functions/FindKmers')
sys.path.append('./Functions/Kmer2Motif')
sys.path.append('./Functions/RanKmers')
sys.path.append('./Functions/TeamKmers')
sys.path.append('./Functions/ViewKmers')
from Functions.FindKmers import FindKmers
from Functions.Kmer2Motif import Kmer2Motif
from Functions.RanKmers import RanKmers
from Functions.TeamKmers import TeamKmers
from Functions.ViewKmers import ViewKmers

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Please provide the module you want to execute. e.g.: python PREDICT.py FindKmers")
        sys.exit(1)

    module = sys.argv[1]

    if module == 'FindKmers':
        parser = argparse.ArgumentParser(description="Run FindKmers Module")
        parser.add_argument('--features', type=str, default='gene', help='Feature type, default is "gene"')
        parser.add_argument('--pThreshold', type=float, default=0.01, help='p-value threshold')
        parser.add_argument('--alg', type=str, default='RandomForest', help='Algorithm name')
        parser.add_argument('--up_stream', type=int, default=1000, help='Upstream length')
        parser.add_argument('--down_stream', type=int, default=500, help='Downstream length')
        args = parser.parse_args(sys.argv[2:])
        FindKmers.main(
            features=args.features,
            pThreshold=args.pThreshold,
            alg=args.alg,
            up_stream=args.up_stream,
            down_stream=args.down_stream
        )

    elif module == 'Kmer2Motif':
        parser = argparse.ArgumentParser(description="Run Kmer2Motif Module")
        parser.add_argument('--TopKmers', type=float, default=0.1, help='Top Kmers proportion')
        parser.add_argument('--KeepTopMotifs', type=float, default=0.1, help='Keep Top Motifs proportion')
        parser.add_argument('--ScoreCutOff', type=float, default=0.9, help='Score cut-off')
        args = parser.parse_args(sys.argv[2:])
        Kmer2Motif.main(
            TopKmers=args.TopKmers,
            KeepTopMotifs=args.KeepTopMotifs,
            ScoreCutOff=args.ScoreCutOff
        )

    elif module == 'RanKmers':
        parser = argparse.ArgumentParser(description="Run RanKmers Module")
        parser.add_argument('--features', type=str, default='gene', help='Feature type')
        parser.add_argument('--alg', type=str, default='RandomForest', help='Algorithm name')
        parser.add_argument('--up_stream', type=int, default=500, help='Upstream length')
        parser.add_argument('--down_stream', type=int, default=1000, help='Downstream length')
        args = parser.parse_args(sys.argv[2:])
        RanKmers.main(
            features=args.features,
            alg=args.alg,
            up_stream=args.up_stream,
            down_stream=args.down_stream
        )

    elif module == 'TeamKmers':
        parser = argparse.ArgumentParser(description="Run TeamKmers Module")
        parser.add_argument('--Motif', type=lambda x: (str(x).lower() == 'true'), default=True, help='Motif flag (True/False)')
        parser.add_argument('--TPCSThreshold', type=float, default=0.75, help='TP cosine similarity threshold(must larger)')
        parser.add_argument('--TPTNCSThreshold', type=float, default=0.3, help='TP&TN cosine similarity threshold(must lower)')
        args = parser.parse_args(sys.argv[2:])
        TeamKmers.main(
            Motif=args.Motif,
            TPCSThreshold=args.TPCSThreshold,
            TPTNCSThreshold=args.TPTNCSThreshold
        )

    elif module == 'ViewKmers':
        parser = argparse.ArgumentParser(description="Run ViewKmers Module")
        parser.add_argument('--features', type=str, default='gene', help='Feature type')
        parser.add_argument('--up_stream', type=int, default=1000, help='Upstream length')
        parser.add_argument('--down_stream', type=int, default=500, help='Downstream length')
        args = parser.parse_args(sys.argv[2:])
        ViewKmers.main(
            features=args.features,
            up_stream=args.up_stream,
            down_stream=args.down_stream
        )

    else:
        print(f"Unknown module: {module}")
