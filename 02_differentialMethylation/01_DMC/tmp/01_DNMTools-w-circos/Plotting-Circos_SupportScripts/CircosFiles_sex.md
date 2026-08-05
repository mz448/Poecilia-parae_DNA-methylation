____
**DATE:**  20250103 
**AUTHOR:** MZF
_____
**GOAL :**
files with configuration parameters necessary for Circos plot
____

# axis.conf
```
axis           = no
axis_color     = dgrey
axis_thickness = 2
axis_spacing   = 0.1
```

# background.conf
```
background                  = transparent
#background_stroke_color     = black
#background_stroke_thickness = 2
```

# bands.conf
```
show_bands            = no
fill_bands            = yes
band_stroke_thickness = 2
band_stroke_color     = white
band_transparency     = 4
```

# Highlight_chromosomes.txt
```bash
awk '{print $1 "\t" 0 "\t" $2}' karyotype.PparFemVer2024.txt > Highlight_chromosomes.txt
```

# ideogram.conf
```
<ideogram>

<spacing>

default = 2u
#break   = 0.01r


axis_break_at_edge = yes
axis_break         = yes
axis_break_style   = 0.1r


<pairwise Parae_23 Parae_01>
spacing = 60r
</pairwise> 


#<break_style 1>
#stroke_color = black
#fill_color   = blue
#thickness    = 0.25r
#stroke_thickness = 2p
#</break>

#<break_style 2>
#stroke_color     = black
#stroke_thickness = 5p
#thickness        = 2r
#</break>

</spacing>

<<include ideogram.position.conf>>
<<include ideogram.label.conf>>
<<include bands.conf>>

</ideogram>
```

# ideogram.label.conf
```
show_label       = yes
label_font       = default
# labels outside circle
#label_radius     = dims(ideogram,radius) + 0.05r
#labels inside circle
#label_radius     = dims(ideogram,radius) - 0.15r
#label inside ideograms
label_radius = (dims(ideogram,radius_inner)+dims(ideogram,radius_outer))/2-24
label_with_tag   = yes
label_size       = 48
label_parallel   = yes
#label_case       = upper
```

# ideogram.position.conf
```
radius           = 0.90r
thickness        = 100p
fill             = yes
fill_color       = black
#stroke_thickness = 2
#stroke_color     = black
```

# config.sex.conf
```
<<include etc/colors_fonts_patterns.conf>>

<<include ideogram.conf>>
#<<include ticks.conf>>

<image>
angle_offset* = -180
<<include etc/image.conf>>
</image>

karyotype = ../karyotype_circos.PparFemVer2024.txt

chromosomes_units           = 1000000
chromosomes_display_default = no

chromosomes_color	    = Parae_01=vlgrey,Parae_02=vlgrey,Parae_03=vlgrey,Parae_04=vlgrey,Parae_05=vlgrey,Parae_06=vlgrey,Parae_07=vlgrey,Parae_08=vlgrey,Parae_09=vlgrey,Parae_10=vlgrey,Parae_11=vlgrey,Parae_12=grey,Parae_13=vlgrey,Parae_14=vlgrey,Parae_15=vlgrey,Parae_16=vlgrey,Parae_17=vlgrey,Parae_18=vlgrey,Parae_19=vlgrey,Parae_20=vlgrey,Parae_21=vlgrey,Parae_22=vlgrey,Parae_23=vlgrey

<highlights>
z = 20
#<highlight>
#file       = Highlight_chromosomes.txt
#r0         = 0.2r
#r1         = 1.01r
#fill_color = white
stroke_color = black
#stroke_thickness = 2u
#</highlight>


<highlight>
#--PARAE REGION 
file	   = Highlight_chromosomes.txt
r0         = 0.71r
r1         = 0.99r
fill_color = white
stroke_color = vlgrey
stroke_thickness = 3u
</highlight>

<highlight>
#--YELLOW REGION
file	   = Highlight_chromosomes.txt
r0         = 0.4r
r1         = 0.7r
fill_color = white
stroke_color = vlgrey
stroke_thickness = 2u
</highlight>

<highlight>
#--IMMACULATA REGION
file	   = Highlight_chromosomes.txt
r0         = 0.2r
r1         = 0.39r
fill_color = white
stroke_color = vlgrey
stroke_thickness = 2u
</highlight>




</highlights>


<plots>

# Make all shared parameters central by including
# them in the outer <plots> block. These values are
# inherited by each <plot> block, where they can
# be further overridden.

type       = histogram
#extend_bin = no
color      = black
#fill_under = yes
thickness*  = 0
min=0
max=0.2

<plot>
#--FEMALE-vs-PARAE HISTOGRAM
file = ../normalized_output.CpGpositions_significant.radadjust_f-vs-p.txt.txt
r0 = 0.71r
r1 = 0.99r
fill_color = purple
</plot>

<plot>
#--FEMALE-vs-YELLOW HISTOGRAM
file = ../normalized_output.CpGpositions_significant.radadjust_f-vs-y.txt.txt
r0 = 0.4r
r1 = 0.7r
fill_color = orange
</plot>

<plot>
#--FEMALE-vs-IMMACULATA HISTOGRAM
file = ../normalized_output.CpGpositions_significant.radadjust_f-vs-i.txt.txt
r0 = 0.2r
r1 = 0.39r
fill_color = grey
</plot>

</plots>
# to explicitly define what is drawn
chromosomes = Parae_02;Parae_03;Parae_01;Parae_09;Parae_16;Parae_05;Parae_07;Parae_17;Parae_13;Parae_10;Parae_04;Parae_06;Parae_18;Parae_12;Parae_11;Parae_15;Parae_14;Parae_08;Parae_22;Parae_19;Parae_21;Parae_20;Parae_23

# to use the logical order ch 1 to 23
chromosomes_order = Parae_01,Parae_02,Parae_03,Parae_04,Parae_05,Parae_06,Parae_07,Parae_08,Parae_09,Parae_10,Parae_11,Parae_12,Parae_13,Parae_14,Parae_15,Parae_16,Parae_17,Parae_18,Parae_19,Parae_20,Parae_21,Parae_22,Parae_23

# To rearange the chromosomes so the 12 is the first in the circle
#chromosomes_order = Parae_12,Parae_1,Parae_2,Parae_3,Parae_4,Parae_5,Parae_6,Parae_7,Parae_8,Parae_9,Parae_10,Parae_11,Parae_13,Parae_14,Parae_15,Parae_16,Parae_17,Parae_18,Parae_19,Parae_20,Parae_21,Parae_22,Parae_23

# To draw  ch12 poping  out!
#chromosomes_radius = Parae_1:0.90r,Parae_2:0.90r,Parae_3:0.90r,Parae_4:0.90r,Parae_5:0.90r,Parae_6:0.90r,Parae_7:0.90r,Parae_8:0.90r,Parae_9:0.90r,Parae_10:0.90r,Parae_11:0.90r,Parae_12:1.00r,Parae_13:0.90r,Parae_14:0.90r,Parae_15:0.90r,Parae_16:0.90r,Parae_17:0.90r,Parae_18:0.90r,Parae_19:0.90r,Parae_20:0.90r,Parae_21:0.90r,Parae_22:0.90r,Parae_23:0.90r

<<include etc/housekeeping.conf>>
```