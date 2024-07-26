	#!/bin/bash
 
	TOKEN=
	CHAT_ID=
	MESSAGE=hello
	
  curl 'https://api.telegram.org/bot$TOKEN/sendMessage?chat_id=$CHAT_ID&text=$MESSAGE'

#### Send Telegram Message ####
curl "https://api.telegram.org/bot6605743120:AAHsqkBpVtDt9vKtJRU7bcsSPLq7RLY3nx0/sendMessage?chat_id=6698838969&text={Finished Static Calculation}"

cp ~/bash/dosplot.py ./
python dosplot.py

curl -s -X POST "https://api.telegram.org/bot6605743120:AAHsqkBpVtDt9vKtJRU7bcsSPLq7RLY3nx0/sendDocument" \
-F "chat_id=6698838969" \
-F "document=@./DOS_Plot.png" \
-F "caption=DOS Plot"
