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
    import subprocess
    import shutil

    # 1. Direct Rclone Google Drive Upload with explicit mimeType for HTML preview
    rclone_bin = shutil.which("rclone") or os.path.expanduser("~/micromamba/envs/cravat_env/bin/rclone")
    if rclone_bin and os.path.exists(rclone_bin):
        cmd = [rclone_bin, "copy", file_path, f"drive:{folder_name}/"]
        if file_path.endswith(".html"):
            cmd.extend(["--header", "Content-Type: text/html"])
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode == 0:
            print(f"[rclone Cloud Sync] Uploaded '{os.path.basename(file_path)}' -> Google Drive: '{folder_name}/{os.path.basename(file_path)}'")
        else:
            print(f"[rclone Warning] {res.stderr.strip()}")

    # 2. Local Google Drive sync directory fallback
    gdrive_local_root = os.path.expanduser("~/Google Drive/My Drive")
    if os.path.exists(gdrive_local_root):
        target_dir = os.path.join(gdrive_local_root, folder_name)
        os.makedirs(target_dir, exist_ok=True)
        dest_file = os.path.join(target_dir, os.path.basename(file_path))
        shutil.copy2(file_path, dest_file)
        print(f"[Local Sync] Copied '{os.path.basename(file_path)}' to '{dest_file}'")

    return True

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
