#!/usr/bin/perl
#
#####################################################################################
# This script works in Grass6.x (hopefully).
#
# The script uses training data from a vector in Grass to create a random forest
# object.  It then creates a multilayer tiff file in Grass, which is tiled and fed
# back into R for the prediction, resulting in the output of class tiles, which are
# read back into Grass and mosaicked.
#
# Mutliscaling is performed on all predictor layers, whose window sizes are stored
# in @scales.  Textural measures are performed at the same scales - the measures
# themselves are stored in %textures.  For each scale, aggregates are calculated
# according to the functions in @aggregates.
#
# Grass rasters to be used in the spectral analysis should be listed in a text file
# called "spectral" in the same directory, and rasters used for the textural 
# analysis should be stored in a file called "textural".  Note with the latter, 
# the raster should be 8-bit integer (CELL type in Grass-speak).  Also, go easy on
# the number of layers - they multiply right up by scales and measures, remember.
#
# You need the R package, with the following packages installed:
# randomForest
# raster
# rgdal
# foreign
#
# It's worth the pain - Random Forest is cool.  By the way, the script tiles the
# multilayer output for the prediction in $tilesize x $tilesize pixel units, and the RF process
# does the prediction in 100,000 pixel batches, both of which measures are designed
# to prevent the process falling over no matter how large the raster files...
#    ...however...
# that's not to say that if your files are too large, the proces won't be in-
# credibly slow, so take care.
#
# What ever your region extents are at the time of invoking the script, those are
# the maximum extents that will be used (and tiled within).  Set your mask carefully,
# it will be maintained.
#
# Cheers
#
# Damien.O'Grady @t astron d0t c0m d0t au
# 11 August 2015
#
# #################################################################################
use strict;
use warnings;

my $tilesize = 500;
my $stage = '';
my $repeat = ($ARGV[0] eq 'repeat') ? 1 : 0;

print "Is the the training (t) stage, or the prediction (p) stage, or in fact the cleanup (c) stage?: ";
while ($stage ne 't' & $stage ne 'p' & $stage ne 'c') {
        $stage = <STDIN>;
        chomp $stage;
}

if ($stage eq 'c') {
	cleanup();
	die "Exiting now\n";
}

my @spred = @{get_layers("spectral")}; # read from file list of layers "spectral"
my @tpred = @{get_layers("textural")}; # read from file list of layers "textural" - say, perhaps just the panchromatic band, or PCA.1 etc
#my %textures = (
#a => 'ASM',
#c => 'Contr',
#k => 'Corr',
#p => 'DE',
#d => 'DV',
#e => 'Entr',
#i => 'IDM',
#s => 'SA',
#w => 'SE',
#x => 'SV',
#v => 'Var'
#); # ALPHABETICAL ORDER!!! hash of r.texture flag => prefix to be used in textural measures
my %textures = (
s => 'SA',
x => 'SV'
); # ALPHABETICAL ORDER!!! hash of r.texture flag => prefix to be used in textural measures
my $training_points = 'training'; # vector file with classes in "class" column
my @scales = (3,7,15,25);
my @tscales = (3);
my $alternate_resolution = undef; # make undef if you don't want to use a larger scale in the textural measures
my @aggregates = ('average','stddev'); # Changing this will break it
my @predictors;
my %codes;

if ($stage eq 't') {
        store_region();
        create_multiscale_layers();
        link_in_dbf_file();
        extract_training_data();
        export_predictors();
	export_codes();
        random_forest();
} elsif ($stage eq 'p') {
        restore_region();
        import_predictors();
        predict();
	cleanup();
}

sub cleanup {
	my $ans = '';
	print "Do you want to delete the temporary files?(y/n): ";
	while ($ans ne 'y' & $ans ne 'n') {
	        $ans = <STDIN>;
	        chomp $ans;
	}
	if ($ans eq 'y') {
		open my $import, "predictors" or die "Couldn't open predictors file.\n";
		while (<$import>) {
			system "g.remove $_ ";
			system "v.db.dropcol training col=$_ ";
		}
	} else {
		print "Nothing deleted.\n";
	}
}

sub restore_region {
	system 'g.region classify ';
}

sub export_codes {
	open my $output, ">codes" or die "Couldn't create codes file.\n";
	foreach my $key (sort keys %codes) {
		print $output "$key:\t$codes{$key}\n";
	}
	close $output;
}

sub import_predictors {
        @predictors=();
        open my $input, "predictors" or die "Couldn't open predictors file.\n";
        while (<$input>) {
                chomp;
                push @predictors, $_;
        }
        close $input;
}

sub export_predictors {
        open my $output, ">predictors" or die "Couldn't create predictors file.\n";
        foreach (@predictors) {
                print $output "$_\n";
        }
        close $output;
}

sub store_region {
        system "g.region save=classify --o ";
}

sub predict {
        print "Carrying out classification now...\n";
        generate_predict_script();
        my $extent_tiles = get_tile_extents();
        make_group();
        my $num_tiles = scalar (keys %{$extent_tiles});
        print "There are $num_tiles tiles.\n";
        my $count = 0;
        foreach my $name (keys %{$extent_tiles}) {
                $count++;
                print "Processing tile $count of $num_tiles...\n";
                my ($n, $s, $e, $w) = @{$extent_tiles->{$name}};
                system "g.region classify; g.region e=$e w=$w n=$n s=$s ";
		my $check = `r.univar -g $predictors[0]`;
		if ($check =~ /^n/) {
	                system "r.out.gdal -c to_classify typ=Float32 out=to_classify.tif ";
	                system "Rscript predict.R";
	                system "r.in.gdal results.tif out=$name --o ";
		} else {
			delete $extent_tiles->{$name};
		}
        }
        mosaic_and_delete_tiles($extent_tiles);
        print "Done\n";
}

sub mosaic_and_delete_tiles {
        my $extent_tiles = shift;
        my @names = keys %$extent_tiles;
        my $string = 'null()';
        foreach my $name (@names) {
                $string =~ s/null\(\)/if(isnull($name),null(),$name)/;
        }
        print "Mosaicking tiles...\n";
	system "g.region classify ";
	system "r.mapcalc 'mosaic=$string' ";
}


sub make_group {
        my $predictor_string = join ",",@predictors;
        system "g.remove group=to_classify ";
        system "i.group to_classify inp=$predictor_string ";
}

sub get_tile_extents {
        my %region = %{get_region()};
        my %output;
        my $num_tile_cols = int($region{cols} / $tilesize) + 1;
        my $num_tile_rows = int($region{rows} / $tilesize) + 1;
        for (my $tile_row = 1; $tile_row <= $num_tile_rows; $tile_row++) {
                my $s = $region{s} + ($tile_row - 1) * $tilesize * $region{nsres};
                my $n = $s + $tilesize * $region{nsres};
                $n = ($n > $region{n}) ? $region{n} : $n;
                for (my $tile_col = 1; $tile_col <= $num_tile_cols; $tile_col++) {
                        my $w = $region{w} + ($tile_col - 1) * $tilesize * $region{ewres};
                        my $e = $w + $tilesize * $region{ewres};
                        $e = ($e > $region{e}) ? $region{e} : $e;

                        $output{"R$tile_row"."C$tile_col"} = [$n, $s, $e, $w];
                }
        }
        return \%output;
}

sub get_region {
        my $result = `g.region -g`;
        my %output = split /[=\n]/, $result;
        return \%output;
}

sub generate_predict_script {
        my $predictor_string = join "','",@predictors;
        $predictor_string = "'$predictor_string'";

        open my $input, ">predict.R" or die "Couldn't create 'predict.R'.\n";
        print $input <<"EOF";
library(raster)
library(rgdal)
library(randomForest)
library(foreign)

rfPredict <- function(rf, imagefile, bandnames, imageout=paste0('p_',imagefile)) {
  o <- raster(imagefile,band=1) #Import the first raster band as prototype
  nbands <- length(bandnames)
  out<-matrix()
  s<-stack(imagefile)
  names(s)<-bandnames

  batches<-ceiling(ncell(s)/100000)
  
  ptm <- proc.time()
  for (i in 1:batches) {
    ulimit <- i * 100000
    llimit <- (ulimit - 100000) + 1
    ulimit <- ifelse(ulimit > ncell(s), ncell(s), ulimit)
    newd<-s[llimit:ulimit]
    newd[is.na(newd)] <- 0
    out[llimit:ulimit]<-predict(rf,newdata=newd,type='prob')[,2] # Now exports probability of DFT
  }
  print(proc.time() - ptm)

  values(o)<-out
  writeRaster(o,imageout,'GTiff',datatype='FLT4S', overwrite=T)
  return(1)
}

load(file = 'rf.RData')

mynames<-c($predictor_string)

fileout='results.tif'
output<-rfPredict(rf,'to_classify.tif',bandnames = mynames, imageout = fileout)

EOF
        close $input;
}


sub random_forest {
        print "Training Random Forest algorithm...\n";
        generate_rfR_script();
        system "Rscript rf.R ";
}


sub extract_training_data {
        print "Extracting training data...\n";
        foreach my $predictor (@predictors) {
                print "Adding $predictor...\n";
                system "v.db.addcol training col='$predictor double' ";
                system "v.what.rast training rast=$predictor col=$predictor ";
        }
}

sub link_in_dbf_file {
        my $result = `g.gisenv -s | sed "s/[;']//g"`;
        my %env = split /[=\n]/,$result;
        my $filename = "$env{GISDBASE}/$env{LOCATION_NAME}/$env{MAPSET}/dbf/training.dbf";
        system "ln -s $filename ";
}

sub generate_rfR_script {
        open my $input, ">rf.R" or die "Couldn't create 'rf.R'.\n";
        print $input <<'EOF';
library(foreign)
library(randomForest)
predictors<-scan('predictors',character())
data<-read.dbf('training.dbf')
classes<-as.factor(data$class)

rf<-randomForest(classes ~ ., data[,predictors], proximity = T, na.action = na.omit)

save(rf, file='rf.RData')

EOF
        close $input;
}

sub create_multiscale_layers {
        print "Creating multiscale layers...\n";
        my $counter = 0;
        foreach my $spectral (@spred) {
		$counter++;
		print "Copying $spectral to B$counter...\n";
		system "g.copy $spectral,B$counter --o ";
		$codes{"B$counter"} = $spectral;
		push @predictors, "B$counter";
                foreach my $scale (@scales) {
			my @twoaggs;
                        foreach my $agg (@aggregates) {
                		$counter++;
                                print "Creating $spectral scale $scale aggregate $agg...\n";
                                system "r.neighbors $spectral out=S$counter met=$agg siz=$scale --o " unless $repeat;
				$codes{"S$counter"} = "$spectral.$scale.$agg";
                                push @predictors, "S$counter";
				push @twoaggs, "S$counter";
                        }
                	$counter++;
                        print "Normalising $spectral scale $scale...\n";
                        system "r.mapcalc 'N$counter = (float($spectral) - $twoaggs[0])/$twoaggs[1] ' " unless $repeat;
			$codes{"N$counter"} = "$spectral.$scale.norm";
                        push @predictors, "N$counter";
                }
        }
        $counter = 0;
        my $flagstring = join ('',keys %textures);
	drop_resolution();
        foreach my $textural (@tpred) {
		my $tex_rescaled = rescale_textural($textural);
                foreach my $tscale (@tscales) {
                	$counter++;
                	$codes{"T$counter"} = "$textural.$tscale";
                        print "creating textures for scale $tscale...\n";
                        system "r.texture -$flagstring $tex_rescaled pre=T$counter --o ";# unless $repeat;
                        add_predictors_to_list("T$counter");
                }
        }
	restore_region();
	truncate_long_names();
}


sub rescale_textural {
	my $textural = shift;
	my $name = "$textural.r";
	if (defined $alternate_resolution) {
		system "r.resamp.stats $textural out=$textural.temp1 met=average --o ";
		system "r.rescale $textural.temp1 out=$textural.temp2 to=0,255 --o";
		system "r.mapcalc '$name = if(isnull($textural.temp2),0,$textural.temp2)' ";
		system "g.remove $textural.temp1,$textural.temp2 -f ";
	} else {
		system "r.rescale $textural out=$textural.temp2 to=0,255 --o";
		system "r.mapcalc '$name = if(isnull($textural.temp2),0,$textural.temp2)' ";
		system "g.remove $textural.temp2 -f ";
	}
	return $name;
}

sub drop_resolution {
	system "g.region res=$alternate_resolution " if defined $alternate_resolution;
}

sub truncate_long_names {
	my @names = `g.mlist rast pat='T*'`;
	chomp @names;
	foreach my $name (@names) {
		if (length($name) > 10) {
			my $newname = substr($name, 0, 10);
			system "g.rename $name,$newname ";
		}
	}
}

sub add_predictors_to_list {
        my $prefix = shift;
        foreach my $measure (sort values %textures) {
                foreach my $angle ((0,45,90,135)) {
                        push @predictors, substr("$prefix"."_$measure"."_$angle",0,10);
                }
        }
}

sub get_layers {
        my $filename = shift;
        open my $input, $filename or die "Couldn't find file '$filename'.\n";
        my @output = (<$input>);
        chomp @output;
	@output = () if ($output[0] eq '') & (scalar(@output)==1);
        return \@output;
}
