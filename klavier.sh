#!/bin/bash
java -Xmx12g -Djava.library.path="./" -jar klavier.jar -inputType FILE -updateActivationsIters 400 -modelNumEMIters 5 -imslpTransModelSer ./data/imslp_indep_transitions.ser -train true -initSpectEnvSer ./data/initSpectEnvFull.ser -spectEnvSer ./data/learnedSpecEnvGould2Itunes.ser -writeTrainedSpect false -filterNoteMinVel 20 -filterNoteMaxIndex 108 -audioFilePath $1
