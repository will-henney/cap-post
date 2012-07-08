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
    "--brightscale", type=float, default=1e7,
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

parser.add_argument(
    "--vcodec", type=str, choices=("wmv2", "x264", "msmpeg4v2"), default="x264",
    help="Video codec to use")

parser.add_argument(
    "--containerformat", type=str, choices=("avi", "mp4"), default="mp4",
    help="Video container format to use")

parser.add_argument(
    "--encoder" type=str, choices=("mencoder", "ffmpeg"), default="ffmpeg",
    help="Video encoding program to use"

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
movie.brightmax = cmd_args.brightscale
movie.bandscales = cmd_args.bandscales
movie.vcodec = cmd_args.vcodec
movie.containerformat = cmd_args.containerformat
movie.encoder = cmd_args.encoder
movie.boxsize = 4.0
movie.camera.set_angles(*cmd_args.orient)

def bmaxNHO(i): 
    """
    Smoothly varying brightness with time so that the movie looks good
    """
    return movie.brightmax/(1.0+(float(i)/10)**2)

def bmaxCPF(i): 
    """
    Smoothly varying brightness with time so that the movie looks good

    This version keeps the R band constant, while G and B bands get brighter with time
    """
    return [movie.brightmax, movie.brightmax/(1.0+(float(i)/70)**2), movie.brightmax/(1.0+(float(i)/50)**2)]



if cmd_args.mode == "tumble":
    movie.camera.set_steps(360.0/movie.nframes, 360.0/movie.nframes)
elif cmd_args.mode == "evo":
    movie.camera.set_steps(0.0, 0.0)
    movie.dtime = 1
    if cmd_args.emshort == "NHO":
        movie.brightmaxfunc = bmaxNHO
    elif cmd_args.emshort == "CPF":
        movie.brightmaxfunc = bmaxCPF
else:
    raise ValueError, "Unknown mode: {}".format(cmd_args.mode)
movie.makemovie()
