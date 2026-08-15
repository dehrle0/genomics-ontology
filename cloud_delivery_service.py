import os
import sys
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders

try:
    from googleapiclient.discovery import build
    from google_auth_oauthlib.flow import InstalledAppFlow
    from google.auth.transport.requests import Request
    from google.oauth2.credentials import Credentials
except ImportError:
    # Handle absence of Google libraries gracefully
    build = None

SCOPES = ['https://www.googleapis.com/auth/drive.file']

def get_gdrive_service():
    if build is None:
        print("[Warning] google-api-python-client is not installed. Skipping Drive connection.")
        return None
    creds = None
    if os.path.exists('token.json'):
        creds = Credentials.from_authorized_user_file('token.json', SCOPES)
    if not creds or not creds.valid:
        if creds and creds.expired and creds.refresh_token:
            creds.refresh(Request())
        else:
            flow = InstalledAppFlow.from_client_secrets_file('config/gdrive_credentials.json', SCOPES)
            creds = flow.run_local_server(port=0)
        with open('token.json', 'w') as token:
            token.write(creds.to_json())
    return build('drive', 'v3', credentials=creds)

def upload_to_drive(file_path: str, folder_name: str = "Ontology"):
    service = get_gdrive_service()
    if not service:
        print("[Warning] Could not initialize Google Drive Service. Save to GDrive skipped.")
        return None

    # Find or create target folder
    query = f"name = '{folder_name}' and mimeType = 'application/vnd.google-apps.folder' and trashed = false"
    results = service.files().list(q=query, fields="files(id)").execute()
    folders = results.get('files', [])

    if not folders:
        # Create folder
        folder_metadata = {
            'name': folder_name,
            'mimeType': 'application/vnd.google-apps.folder'
        }
        folder = service.files().create(body=folder_metadata, fields='id').execute()
        folder_id = folder.get('id')
        print(f"Created new GDrive folder: '{folder_name}' (ID: {folder_id})")
    else:
        folder_id = folders[0]['id']
        print(f"Found existing GDrive folder: '{folder_name}' (ID: {folder_id})")

    # Upload file
    file_name = os.path.basename(file_path)
    file_metadata = {
        'name': file_name,
        'parents': [folder_id]
    }
    
    # Simple media upload
    from googleapiclient.http import MediaFileUpload
    media = MediaString = MediaFileUpload(file_path, resumable=True)
    uploaded_file = service.files().create(body=file_metadata, media_body=media, fields='id').execute()
    
    file_id = uploaded_file.get('id')
    print(f"Successfully uploaded '{file_name}' to GDrive folder '{folder_name}' (File ID: {file_id})")
    return file_id

def send_email_notification(patient_id: str, report_filepath: str, recipient_email: str):
    sender_email = os.environ.get("SENDER_EMAIL", "pipeline@genomics-ontology.org")
    sender_password = os.environ.get("SENDER_PASSWORD", "")
    smtp_server = os.environ.get("SMTP_SERVER", "smtp.gmail.com")
    smtp_port = int(os.environ.get("SMTP_PORT", "587"))

    if not recipient_email:
        print("[Warning] No recipient email configured. Skpping email dispatch.")
        return

    msg = MIMEMultipart()
    msg['From'] = sender_email
    msg['To'] = recipient_email
    msg['Subject'] = f"Genomic Report Completed: {patient_id}"

    body = f"Hello,\n\nThe whole-genome sequencing report for patient {patient_id} has been successfully completed and uploaded to Google Drive.\n\nBest regards,\nGenomics Pipeline Automation"
    msg.attach(MIMEText(body, 'plain'))

    # Attach HTML report file
    try:
        with open(report_filepath, "rb") as attachment:
            part = MIMEBase("application", "octet-stream")
            part.set_payload(attachment.read())
            encoders.encode_base64(part)
            part.add_header(
                "Content-Disposition",
                f"attachment; filename= {os.path.basename(report_filepath)}",
            )
            msg.attach(part)
    except Exception as e:
        print(f"[Error] Failed to attach report file: {e}")

    # Send mail
    try:
        server = smtplib.SMTP(smtp_server, smtp_port)
        server.starttls()
        if sender_password:
            server.login(sender_email, sender_password)
        text = msg.as_string()
        server.sendmail(sender_email, recipient_email, text)
        server.quit()
        print(f"Successfully sent completion email to {recipient_email}")
    except Exception as e:
        print(f"[Warning] SMTP Connection failed: {e}. Email notification could not be dispatched.")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 cloud_delivery_service.py <file_path> <patient_id>")
        sys.exit(1)
    file_path = sys.argv[1]
    patient_id = sys.argv[2]
    upload_to_drive(file_path)
