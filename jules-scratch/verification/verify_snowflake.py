import os
import time
from playwright.sync_api import sync_playwright

def run():
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        page = browser.new_page()

        # Go to the local HTML file via HTTP server
        page.goto('http://localhost:8000/index.html')

        # Wait for the scene to be rendered
        page.wait_for_selector('canvas')

        # Give it a moment to run the simulation
        time.sleep(5)

        # Take a screenshot
        page.screenshot(path='jules-scratch/verification/verification.png')

        browser.close()

if __name__ == '__main__':
    run()
