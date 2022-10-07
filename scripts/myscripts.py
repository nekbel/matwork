##########################
# python mail sent
##########################

import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText

mail_content = "ΞΥΠΝΑ ΕΛΛΗΝΑ Ο ΓΙΟΣ ΣΟΥ ΤΟΝ ΒΡΟΝΤΑΕΙ ΜΕ ΚΙΝΕΖΙΚΑ ΜΙΚΙ ΜΑΟ"

#The mail addresses and password
sender_address = 'zorznekbel100@gmail.com'
sender_pass = 'xbtcvyiiinzbjypz'
receiver_address = 'nek.belmezos@gmail.com'
#Setup the MIME
message = MIMEMultipart()
message['From'] = sender_address
message['To'] = receiver_address
message['Subject'] = 'Script done running'   #The subject line
#The body and the attachments for the mail
message.attach(MIMEText(mail_content, 'plain'))
#Create SMTP session for sending the mail
session = smtplib.SMTP('smtp.gmail.com', 587) #use gmail with port
session.starttls() #enable security
session.login(sender_address, sender_pass) #login with mail_id and password
text = message.as_string()
session.sendmail(sender_address, receiver_address, text)
session.quit()
print('Mail Sent')

##########################
# arguments to .py script
##########################
import sys

def hello(a,b):
    print "hello and that's your sum:", a + b

if __name__ == "__main__":
    a = int(sys.argv[1])
    b = int(sys.argv[2])
    hello(a, b)