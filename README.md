# Optimal Tuning Systems

This is an implementation of the Optimal Tuning System algorithm outlined in
[Polansky et al,
2008](http://eamusic.dartmouth.edu/~larry/published_articles/owt_pnm.pdf).

## Build

    make # builds owt
    make key_weights # builds the key_weights search program

## Specify Criteria

    # criteria.owt
    title: HAHAHA 
    num_pitches: 12
    ideal_intervals: 100 204 267 386 498 600 702 800 900 969 1100
    interval_weights: 0.1 0.1 0.1 100.0 0.1 0.1 0.1 0.1 0.1 0.1 0.1
    key_weights: 100 0.1 0.1 0.1 0.1 100 0.1 100 0.1 0.1 0.1 0.1
    repeat_factor: 1200.0

## Run

    # finds a temperament based on the criteria
    ./owt -i criteria.owt 
    # finds key weights based on interval weights and a target temperament
    ./key_weights -i werck_iii.owt 

