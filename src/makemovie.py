from movierot import Movie 
import argparse
import sys

parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description="""
Generate color movies from emissivity maps (evolution and tumble options)
    """)

parser.add_argument(
    "runid", type=str,
    help="ID of simulation run to use (e.g., 04052012_4")
    
parser.add_argument(
    "--mode", type=str, choices=("tumble", "evo"), default="tumble",
    help="Type of movie to make")

parser.add_argument(
    "--time", type=int, default=1,
    help="Time to use (for tumble) or start time (for evo)")

parser.add_argument(
    "--brightscale", type=float, default=2.e8,
    help="Scale for maximum brightness in the middle channel")

parser.add_argument(
    "--bandscales", type=float, nargs=3, default=[0.28, 1.0, 0.2],
    help="Relative adjustments to scales in the three RGB channels")

parser.add_argument(
    "--emtypes", type=str, nargs=3, default=Movie.emtypes,
    help="Emission types to use for the three RGB channels")

parser.add_argument(
    "--emshort", type=str, default="NHO",
    help="Short identifier for the emission types")

parser.add_argument(
    "--orient", type=float, nargs=2, metavar=("THETA", "PHI"), default=[0.0, 0.0],
    help="Initial orientation of cube")

parser.add_argument(
    "--frames", type=int, default=72,
    help="Number of frames in movie")


cmd_args = parser.parse_args()
       
movieid = "{}-t{:04}-th{}-ph{}-n{}".format(cmd_args.mode, cmd_args.time, 
                                           cmd_args.orient[0], cmd_args.orient[1], 
                                           cmd_args.frames)
movie = Movie(
    runid=cmd_args.runid, 
    movieid=movieid,
    time=cmd_args.time,
    nframes=cmd_args.frames,
    verbose=1, force=0,
    )
movie.imageprefix = "rgb-{}-{}-{}".format(cmd_args.emshort, movie.runid, movie.movieid)
movie.emtypes = cmd_args.emtypes
movie.bandscales = cmd_args.bandscales
movie.boxsize = 4.0
movie.camera.set_angles(*cmd_args.orient)
movie.camera.set_steps(360.0/movie.nframes, 360.0/movie.nframes)
movie.makemovie()
