import sys
import os
import argparse

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Please provide the module you want to execute. e.g.: python PREDICT.py KmerImpact")
        sys.exit(1)

    module = sys.argv[1]

    if module == 'KmerImpact':
        sys.path.append('./Functions/KmerImpact')
        from Functions.KmerImpact import KmerImpact
        parser = argparse.ArgumentParser(description="Run KmerImpact Module")
        ## Required
        parser.add_argument('--gff', type=str, required=True, help='Path of your gff file')
        parser.add_argument('--genome', type=str, required=True, help='Path of your genome file')
        parser.add_argument('--tp', type=str, required=True, help='Path of your tp genes file')
        parser.add_argument('--tn', type=str, required=True, help='Path of your tn genes file')

        ## Customized
        parser.add_argument('--features', type=str, default='gene', help='Feature type, default is "gene"')
        parser.add_argument('--pThreshold', type=float, default=0.01, help='p-value threshold')
        parser.add_argument('--alg', type=str, default='RandomForest', help='Algorithm name')
        parser.add_argument('--up_stream', type=int, default=1000, help='Upstream length')
        parser.add_argument('--down_stream', type=int, default=500, help='Downstream length')
        args = parser.parse_args(sys.argv[2:])
        KmerImpact.main(
            gff=args.gff,
            genome=args.genome,
            tp=args.tp,
            tn=args.tn,
            features=args.features,
            pThreshold=args.pThreshold,
            alg=args.alg,
            up_stream=args.up_stream,
            down_stream=args.down_stream
        )

    elif module == 'Kmer2Motif':
        sys.path.append('./Functions/Kmer2Motif')
        from Functions.Kmer2Motif import Kmer2Motif
        parser = argparse.ArgumentParser(description="Run Kmer2Motif Module")
        ## Required
        parser.add_argument('--kmer', type=str, required=True, help='Path of your kmer list file')
        parser.add_argument('--motif', type=str, required=True, help='Path of your motif list file')

        ## Customized
        parser.add_argument('--TopKmers', type=float, default=0.1, help='Top Kmers proportion')
        parser.add_argument('--KeepTopMotifs', type=float, default=0.1, help='Keep Top Motifs proportion')
        parser.add_argument('--ScoreCutOff', type=float, default=0.9, help='Score cut-off')
        args = parser.parse_args(sys.argv[2:])
        Kmer2Motif.main(
            Kmer=args.kmer,
            Motif=args.motif,
            TopKmers=args.TopKmers,
            KeepTopMotifs=args.KeepTopMotifs,
            ScoreCutOff=args.ScoreCutOff
        )

    elif module == 'RanKmers':
        sys.path.append('./Functions/RanKmers')
        from Functions.RanKmers import RanKmers
        parser = argparse.ArgumentParser(description="Run RanKmers Module")
        ## Required
        parser.add_argument('--gff', type=str, required=True, help='Path of your gff file')
        parser.add_argument('--genome', type=str, required=True, help='Path of your genome file')
        parser.add_argument('--tp', type=str, required=True, help='Path of your tp genes file')
        parser.add_argument('--tn', type=str, required=True, help='Path of your tn genes file')
        parser.add_argument('--kmer', type=str, required=True, help='Path of your kmer list file')

        ## Customized
        parser.add_argument('--features', type=str, default='gene', help='Feature type')
        parser.add_argument('--alg', type=str, default='RandomForest', help='Algorithm name')
        parser.add_argument('--up_stream', type=int, default=500, help='Upstream length')
        parser.add_argument('--down_stream', type=int, default=1000, help='Downstream length')
        args = parser.parse_args(sys.argv[2:])
        RanKmers.main(
            gff=args.gff,
            genome=args.genome,
            tp=args.tp,
            tn=args.tn,
            Kmer=args.kmer,
            features=args.features,
            alg=args.alg,
            up_stream=args.up_stream,
            down_stream=args.down_stream
        )

    elif module == 'TeamKmers':
        sys.path.append('./Functions/TeamKmers')
        from Functions.TeamKmers import TeamKmers
        parser = argparse.ArgumentParser(description="Run TeamKmers Module")
        ## Required
        parser.add_argument('--tp', type=str, required=True, help='Path of your (tp)geneXkmer file')
        parser.add_argument('--tn', type=str, required=True, help='Path of your (tn)geneXkmer file')

        ## Customized
        parser.add_argument('--motif', type=str, default="", help='Path of your motif list file')
        parser.add_argument('--TPCSThreshold', type=float, default=0.75, help='TP cosine similarity threshold(must larger)')
        parser.add_argument('--TPTNCSThreshold', type=float, default=0.3, help='TP&TN cosine similarity threshold(must lower)')
        args = parser.parse_args(sys.argv[2:])
        TeamKmers.main(
            tp=args.tp,
            tn=args.tn,
            Motif=args.motif,
            TPCSThreshold=args.TPCSThreshold,
            TPTNCSThreshold=args.TPTNCSThreshold
        )

    elif module == 'ViewKmers':
        sys.path.append('./Functions/ViewKmers')
        from Functions.ViewKmers import ViewKmers
        parser = argparse.ArgumentParser(description="Run ViewKmers Module")
        ## Required
        parser.add_argument('--gff', type=str, required=True, help='Path of your gff file')
        parser.add_argument('--genome', type=str, required=True, help='Path of your genome file')
        parser.add_argument('--gene', type=str, required=True, help='Path of your gene list file')
        parser.add_argument('--kmer', type=str, required=True, help='Path of your kmer list file')

        ## Customized
        parser.add_argument('--motif', type=str, default="", help='Path of your motif list file')
        parser.add_argument('--features', type=str, default='gene', help='Feature type')
        parser.add_argument('--up_stream', type=int, default=1000, help='Upstream length')
        parser.add_argument('--down_stream', type=int, default=500, help='Downstream length')
        args = parser.parse_args(sys.argv[2:])
        ViewKmers.main(
            gff=args.gff,
            genome=args.genome,
            gene=args.gene,
            Kmer=args.kmer,
            Motif=args.motif,
            features=args.features,
            up_stream=args.up_stream,
            down_stream=args.down_stream
        )

    else:
        print(f"Unknown module: {module}")
