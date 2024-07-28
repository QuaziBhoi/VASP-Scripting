#!/bin/bash
 
TOKEN=6605743120:AAHsqkBpVtDt9vKtJRU7bcsSPLq7RLY3nx0
CHAT_ID=6698838969
MESSAGE="{Finished Calculation}"
 
curl "https://api.telegram.org/bot$TOKEN/sendMessage?chat_id=$CHAT_ID&text=$MESSAGE"
echo 

#### Send Telegram Message ####
#curl "https://api.telegram.org/bot6605743120:AAHsqkBpVtDt9vKtJRU7bcsSPLq7RLY3nx0/sendMessage?chat_id=6698838969&text={Finished Static Calculation}"
