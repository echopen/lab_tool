# echOpen lab tool
echOpen GUI software lab tool to test and develop signal processing algorithms. 

## To start
### info
- Software tool [design brief](https://github.com/echopen/lab_tool/blob/master/DB_laboratory_tool_v2_en.pdf)
- Joining the team on echOpen slack [lab tool dedicated channel](https://echopen.slack.com/messages/CCGEF6CQY/)

### Contact (lead dev)
- Jérôme (@GG23800)
- Clément (@clecoued)

## Architecture
The architecture of this software is described in figure below, it is composed of 7 principal modules:

- Connexion management: module that checks the connexion between the probe and the display software.
- Settings management: module that manages the settings between the probe and the display module.
- Feedback: module that manages the feedback from the probe.
- Data streaming reception: module that receives the data from the probe sends it the the signal processing module.
- Signal processing: module that applies the signal processing to each measurement line and sends it to image processing module.
- Image processing: module that determines the image.
- Image display: main window of the software where the image is displayed.
- The black lines represent the flux of the data from the probe to the displayed image. The gray lines represent the flux of settings and some informations such as feedback from the probe.

![architecture scheme](img/archi.png)

## Installing (coming soon...) 
