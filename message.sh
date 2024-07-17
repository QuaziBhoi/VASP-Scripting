	#!/bin/bash
 
	TOKEN=6605743120:AAHsqkBpVtDt9vKtJRU7bcsSPLq7RLY3nx0
	CHAT_ID=6698838969
	MESSAGE=hello
	
  curl 'https://api.telegram.org/bot$TOKEN/sendMessage?chat_id=$CHAT_ID&text=$MESSAGE'
