globals
[ 
  test
  HSL
  EnergyG
]

breed [ cell ]

patches-own
[
]

cell-own
[
  Energy
  LuxR
  LuxI
  Complex
  flag
  duplic_time       
]

;;;;;;;;;;;
;; Setup ;;
;;;;;;;;;;;

to setup
  ca
  
  if Seed?
    [ random-seed 47823 ]
  
  repeat 100 [ cel ]
  
  set HSL 100
  
  set EnergyG Global
  
  ask turtles [ graphic ]
  
end

;;;;;;;;;;;;;;;;;;;;;;;
;; globals procedure ;;
;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;
;; patch procedure ;;
;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;
;; turtles procedure ;;
;;;;;;;;;;;;;;;;;;;;;;;

to cel
  set-default-shape cell "circle 2" 
  crt 1 [ setxy random-pxcor random-pycor 
          set color white
          set breed cell
          set Energy 1000
          set duplic_time 30 + random 30 
          set LuxR 1
          set LuxI 1
          set Complex 0
          set flag 0
          set size 2  
        ]    
end

to move
  let neighbor [ list pxcor pycor ] of neighbors
  let value one-of neighbor
  if ( not any? turtles-on patch ( item 0 value ) ( item 1 value ) )
    [ move-to patch ( item 0 value ) ( item 1 value ) ]  
end

to lumini
  set color scale-color white Complex 0 10
end

;;;;;;;;;;;;;;
;; Run Time ;;
;;;;;;;;;;;;;;

to go
  tick            
  ask cell 
    [ move ]
  
  ;; Evaporation Rate     
  set HSL HSL * ( 1 - 0.10 )
  
  ask cell
    [  ifelse ( Energy < ( ( DeltaE * LuxI-rate ) + 1 ) )
         [ if ( count cell > 100 )
             [ die ]
         ]    
         [ set Energy Energy - 1                           ;; Perda de energia no passo (manutencao)
           set Energy Energy - ( DeltaE * LuxI-rate )      ;; Perda de energia pelo efeito da luciferase
                 
           if ( EnergyG > ( ( DeltaE * LuxI-rate ) + 1 ) )
             [ set EnergyG EnergyG - ( ( DeltaE * LuxI-rate ) + 1 )
               set Energy Energy + ( ( DeltaE * LuxI-rate ) + 1 )
             ]      
           ] 
     ]
     
  ifelse ( precision ( HSL / 40401 ) 2 <= 0.20 )
    [ ;; Low Density
    
      ask cell [ set Complex 0 ]
      
      set EnergyG EnergyG + 10000
      
      set test 0
      
      ask cell  
        [ if ( ticks mod ( duplic_time ) = 0 )
            [ ifelse ( flag < 3 ) 
                [ ifelse ( Energy >= ( count cell ) ) ;;( DeltaE * LuxI ) ) 
                    [ set Energy Energy - 100
                      hatch-cell 1 ;; Replication cellular
                        [ setxy random-pxcor random-pycor 
                          set color white
                          set Energy 100
                          set duplic_time 30 + random 30
                          set LuxR 1
                          set LuxI 1
                          set Complex 0
                          set flag 0
                          set size 2   
                        ]
                      set flag 0   
                    ]
                    [ ifelse ( EnergyG >= ( count cell ) ) ;;( DeltaE * LuxI ) ) 
                        [ set EnergyG EnergyG - 100
                          hatch-cell 1 ;; Replication cellular
                            [ setxy random-pxcor random-pycor 
                              set color white
                              set Energy 100
                              set duplic_time 30 + random 30
                              set LuxR 1
                              set LuxI 1
                              set Complex 0
                              set flag 0
                              set size 2   
                            ]
                          set flag 0        
                        ]
                        [ set flag flag + 1 ]
                  ] 
               ]
               [ die ] 
           ]
         set HSL HSL + HSL-rate
         set LuxR LuxR + LuxR-rate
         set LuxI LuxI + LuxI-rate        
        ]
    ]                
    [ ;; High Density
      ask cell 
        [ lumini
          set test 1 
          if ( ticks mod ( duplic_time ) = 0 )
            [ ifelse ( flag < 3 ) 
                [ ifelse ( Energy >= ( count cell ) ) ;;( DeltaE * LuxI ) ) 
                    [ set Energy Energy - 100
                      hatch-cell 1 ;; Replication cellular
                        [ setxy random-pxcor random-pycor 
                          set color white
                          set Energy 100
                          set duplic_time 30 + random 30
                          set LuxR 1
                          set LuxI 1
                          set Complex 0
                          set flag 0
                          set size 2   
                        ]
                       set flag 0         
                    ] 
                    [ ifelse ( EnergyG >= ( count cell ) ) ;;( DeltaE * LuxI ) ) 
                        [ set EnergyG EnergyG - 100 ;; Retira da energia do sistema
                          hatch-cell 1 ;; Replication cellular
                            [ setxy random-pxcor random-pycor 
                              set color white
                              set Energy 100
                              set duplic_time 30 + random 30
                              set LuxR 1
                              set LuxI 1
                              set Complex 0
                              set flag 0
                              set size 2   
                            ]
                          set flag 0       
                        ]
                        [ set flag flag + 1 ]
                   ] 
               ]
               [ die ]  
                
            ]  
          set HSL HSL + ( 2 * HSL-rate )
          set LuxI LuxI + ( 2 * LuxI-rate ) 
          set LuxR LuxR + ( LuxR-rate / 2 ) 
          if ( ( HSL > 0 ) and ( LuxR > 0 ) )
            [ set Complex Complex + 1 ;; Luminescencia
              set LuxR LuxR - 1
              set HSL HSL - 1
            ]  
        ]   
    ]
    
      
    ;;set-current-plot-pen "HSL"
    ;;plot sum [ HSL ] of patches  
    ;;set-current-plot-pen "LuxR"
    ;;plot sum [ LuxR ] of cell
    set-current-plot-pen "LuxI"
    plot sum [ LuxI ] of cell
    set-current-plot-pen "Complex"
    plot sum [ Complex ] of cell
    
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Plotting procedures ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;

to graphic
  ;; plot the gaphics
  set-current-plot "Dados"
  ;; we don't want the "auto-plot" feature to cause the
  ;; plot's x range to grow when we draw the axis.  so
  ;; first we turn auto-plot off temporarily   
  set-plot-x-range 0 1
  set-plot-y-range 0 1 
  auto-plot-off
  auto-plot-on
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; I/O procedures      ;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;

to images
  let cont 0
  if cont < 10
   [ export-view ( word "0000000" cont ".png" ) ]
  if cont > 10 and cont < 100
   [ export-view ( word "000000" cont ".png" ) ]
  if cont > 100 and cont < 1000
   [ export-view ( word "00000" cont ".png" ) ]  
  set cont cont + 1
end

;
; (C) 2004 Uri Wilensky. From where if it originated the present part of code here.
;
; (C) 2007 Luis C. da Costa
; This code may be freely copied, distributed, altered, 
; or otherwise used by anyone for any legal purpose.
; 
; This model was created as part of the project:
; xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx.
; The project gratefully acknowledges the support of the
; National Laboratory for Scientific Computing - NLSC 
; Computer Science Coordination - CSC
; Petropolis - Rio de Janeiro - Brazil.
;
; THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
; "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
; LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
; A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
; OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
; SPECIAL, EXEMPLARY, OR CONSEquimuimUENTIAL DAMAGES (INCLUDING, BUT NOT
; LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
; DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
; THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
; (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
; OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.    
@#$#@#$#@
GRAPHICS-WINDOW
544
10
1157
644
100
100
3.0
1
10
1
1
1
0
1
1
1
-100
100
-100
100
0
0
1
ticks

CC-WINDOW
5
713
1166
808
Command Center
0

BUTTON
40
21
167
54
NIL
Setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL

BUTTON
200
22
327
55
Execute
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL

SWITCH
30
151
169
184
Seed?
Seed?
1
1
-1000

MONITOR
409
447
476
492
Hours
ceiling ( ticks / 30 )
17
1
11

PLOT
21
207
433
425
Dados
NIL
NIL
0.0
10.0
0.0
10.0
false
true
PENS
"LuxI" 1.0 0 -2674135 true
"Complex" 1.0 0 -16777216 true

MONITOR
18
449
171
494
Cellular Energy 
mean [ Energy ] of cell
17
1
11

MONITOR
223
510
374
555
LuxR
sum [ LuxR ] of cell
17
1
11

MONITOR
223
448
375
493
LuxI
sum [ LuxI ] of cell
17
1
11

MONITOR
20
512
172
557
HSL
HSL
17
1
11

MONITOR
20
581
172
626
General Energy
EnergyG
17
1
11

MONITOR
222
579
377
624
Cellula
count cell
17
1
11

MONITOR
21
653
300
698
NIL
precision ( HSL / 40401) 2
17
1
11

MONITOR
334
650
391
699
NIL
test
17
1
12

SLIDER
346
64
481
97
LuxI-rate
LuxI-rate
0
1.0
1
0.1
1
NIL
HORIZONTAL

SLIDER
345
106
482
139
LuxR-rate
LuxR-rate
0
1.0
1
0.1
1
NIL
HORIZONTAL

SLIDER
345
23
483
56
HSL-rate
HSL-rate
0
10.0
1
1.0
1
NIL
HORIZONTAL

SLIDER
189
151
327
184
DeltaE
DeltaE
0
10
5
1
1
NIL
HORIZONTAL

SLIDER
346
151
484
184
Global
Global
0
100000
1000
100
1
NIL
HORIZONTAL

@#$#@#$#@
WHAT IS IT?
-----------

The discovery that bacteria are able to communicate with each other
changed our general perception of many single, simple organisms
inhabiting our world. Instead of language, bacteria use signalling
molecules which are released into the environment. As well as
releasing the signalling molecules, bacteria are also able to
measure the number (concentration) of the molecules within a
population. Nowadays we use the term 'Quorum Sensing' (QS) to
describe the phenomenon whereby the accumulation of signalling
molecules enables a single cell to sense the number of bacteria
(cell density). In natural environments, there are many types of
bacteria with a variety of signalling molecules. As they employ
different languages they cannot necessarily talk to all other
bacteria. Presenting the model of cellular automata proposed to
describe the main mechanisms of Quorum Sensing where Vibrio
fischeri, and its model using the concept of Multi-Agents System.


HOW TO USE IT
-------------

The basic controls for the model are: 
SETUP - Creating the environment and distribution of Vibrio Fisheri
GO - Run the model 
SLIDERS - Parameters of the experiment
CHOOSER - Mode Display
SWITCH - Selecting a seed


EXTENDING THE MODEL
-------------------

Show potential of modeling multi-agents in models that are repeated in biology.


CREDITS AND REFERENCES
----------------------

J. Liu. Autonomous agents and multi-agent systems: explorations in learning, self- organization and adaptive computation. World Scientiﬁc, 2001.

E. Greenberg. Quorum sensing in gram-negative bacteria. ASM News, 63:371–377, 1997.

P. C. Chen. A computational model of a class of gene networks with positive and negative controls. Technical report, Faculty of Engineering - National University of Singapore, www.elsevier.com/locate/biosystems, 2003.

To refer to this model in academic publications, please use:  Wilensky, U. (1998).  NetLogo Voting model.  http://ccl.northwestern.edu/netlogo/models/Voting.  Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

In other publications, please use:  Copyright 1998 Uri Wilensky.  All rights reserved.  See http://ccl.northwestern.edu/netlogo/models/Voting for terms of use.

@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 4.0.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
