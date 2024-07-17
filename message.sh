	#!/bin/bash
 
	TOKEN=
	CHAT_ID=
	MESSAGE=hello
	
  curl 'https://api.telegram.org/bot$TOKEN/sendMessage?chat_id=$CHAT_ID&text=$MESSAGE'
