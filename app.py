import os
from flask import Flask, request, jsonify, render_template, send_file
import time
import re
import requests
from functools import lru_cache
from typing import Dict, Any, Optional, Tuple, List
from Bio.Seq import Seq
import pandas as pd
import io

# --- Template Setup ---
# This section will automatically create the necessary HTML files in a 'templates' folder.

def setup_templates():
    """Creates the templates directory and HTML files if they don't exist."""
    if not os.path.exists('templates'):
        os.makedirs('templates')

    # Base template with navigation and footer
    base_html = """
<!doctype html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <title>{{ title }} - AVEC</title>
    <script src="https://cdn.jsdelivr.net/npm/igv@2.15.5/dist/igv.min.js"></script>
    <style>
        body{font-family:system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,"Helvetica Neue",Arial,sans-serif;margin:0;background-color:#f8f9fa;color:#333}
        .container{max-width:900px;margin:auto;background:var(--card);padding:2em;border-radius:8px;box-shadow:0 4px 8px rgba(0,0,0,.1); margin-top: 2em; margin-bottom: 2em;}
        h3,h4,h5{color:#005a9c} a{color:#005a9c; text-decoration: none;} a:hover{text-decoration: underline;}
        nav {background-color: #333; padding: 1em; text-align: center;}
        /* THIS IS THE ONLY CHANGE IN THIS BLOCK: Added API link */
        nav a {color: white; margin: 0 15px; text-decoration: none; font-weight: bold;}
        .disclaimer{font-size:.8em;color:#6c757d;margin-top:2em;text-align:center; padding-top: 1em; border-top: 1px solid #eee;}
        form{display:grid;grid-template-columns:1fr;gap:1em;align-items:center;margin-bottom:2em}
        input,button{padding:.75em;border-radius:4px;border:1px solid #ccc;font-size:1em;width:100%;box-sizing:border-box}
        button{background-color:#007bff;color:white;font-weight:bold;cursor:pointer;border:none;width:auto;justify-self:end;padding:.75em 1.5em}
        button:hover{background-color:#0056b3}
        #loader{display:none;text-align:center;padding:1em;font-size:1.2em}
        #results{display:none}.result-header{padding:1em;color:white;border-radius:4px 4px 0 0;margin-bottom:0;}.result-header h4{margin:0;font-size:1.5em;text-align:center;color:inherit;}
        .result-block{border: 1px solid #ddd; border-top: none; padding: 1em; margin-bottom: 1.5em; border-radius: 0 0 4px 4px;}.result-block h5{margin-bottom:.5em;color:#495057;border-bottom:1px solid #eee;padding-bottom:.3em; margin-top: 0;}
        .checklist{list-style-type:none;padding:0}.checklist li{margin-bottom:.5em}.pass::before{content:'\2713';color:#28a745;margin-right:10px;font-weight:bold}.fail::before{content:'\2717';color:#dc3545;margin-right:10px;font-weight:bold}
        .eligible{background-color:#28a745}.likely-eligible{background-color:#17a2b8}.unlikely-eligible{background-color:#ffc107;color:#333 !important;}.not-eligible{background-color:#dc3545}.unable-to-assess{background-color:#6c757d}.error{background-color:#dc3545}.note{font-style:italic;color:#555;font-size:.9em}
        .strategy-block {margin-bottom: 2em;}
        .summary-block { background-color: #e9ecef; padding: 1em; border-radius: 4px; margin-bottom: 2em; }
        #igv-container { border: 1px solid #ddd; margin-top: 2em; margin-bottom: 1.5em; }
        /* Style for code blocks in API docs */
        pre { background-color: #f0f0f0; border: 1px solid #ddd; border-radius: 4px; padding: 1em; white-space: pre-wrap; word-wrap: break-word; font-family: 'Courier New', Courier, monospace; }
        /* PubMed results list */
        .pubmed-list { list-style: none; margin: .5em 0 0 0; padding: 0; display: grid; gap: .75em; }
        .pubmed-card { border: 1px solid var(--border); border-radius: 8px; padding: .75em .9em; background: #fff; box-shadow: 0 1px 2px rgba(0,0,0,.03); }
        .pubmed-title { font-weight: 600; margin: 0 0 .25em; color: #0f172a; }
        .pubmed-meta { color: #555; font-size: 0.9em; }
        .pubmed-card-footer { margin-top: .6em; text-align: right; }
        .pubmed-article-link { display:inline-block; padding: .35em .6em; border: 1px solid var(--border); border-radius: 6px; text-decoration: none; color: #0c63bd; background:#f8fafc }
        .pubmed-article-link:hover { background:#eef2f7 }
        .pubmed-footer { margin-top: .75em; text-align: center; }
        .pubmed-footer a { text-decoration: none; padding:.4em .75em; border-radius:6px; border:1px solid var(--border); color:#0c63bd; background:#fff }
        .pubmed-footer a:hover { background:#f8fafc }

        /* Beauty overrides and enhancements */
        :root{--primary:#0d6efd;--primary-dark:#0b5ed7;--text:#333;--muted:#6c757d;--bg:#f7f8fb;--card:#ffffff;--border:#e6e9ef;--ok:#28a745;--info:#17a2b8;--warn:#ffc107;--bad:#dc3545;--neutral:#6c757d}
        body{background:linear-gradient(180deg,#f8fafc 0%,#eef2f7 100%)}
        .container{max-width:980px;border-radius:12px;box-shadow:0 6px 18px rgba(0,0,0,.08)}
        nav{background:#222;box-shadow:0 2px 6px rgba(0,0,0,.15);padding:.85em 1em}
        nav{position:relative;display:flex;align-items:center;justify-content:center}
        nav .links{display:flex;align-items:center;gap:12px}
        nav .links a{opacity:.9;font-weight:600;margin:0 12px;color:#f3f3f3;text-decoration:none}
        nav .links a:hover{opacity:1}
        .theme-toggle{position:absolute;right:1rem;top:50%;transform:translateY(-50%);background:transparent;border:1px solid var(--border);color:#f3f3f3;border-radius:999px;padding:.35em .7em;cursor:pointer}
        .theme-toggle:hover{background:rgba(255,255,255,.06)}
        .disclaimer{font-size:.85em;border-top:1px solid var(--border)}
        #assessment-form{grid-template-columns:120px 1fr auto}
        #assessment-form label{margin:0;font-weight:600;color:#495057}
        input,button{border-radius:6px;border:1px solid #cfd4dc}
        button{background:var(--primary);box-shadow:0 2px 8px rgba(13,110,253,.25)}
        button:hover{background:var(--primary-dark)}
        .result-header{border-radius:10px 10px 0 0}
        .result-header h4{font-size:1.2em;letter-spacing:.2px}
        .accordion-toggle .chevron{display:inline-block;margin-left:.4em;transition:transform .2s ease}
        .accordion-toggle.open .chevron{transform:rotate(180deg)}
        /* Rounded corners when collapsed; squared bottom when open */
        .accordion-toggle{border-radius:10px}
        .accordion-toggle.open{border-radius:10px 10px 0 0}
        .result-block{border:1px solid var(--border);border-top:none;border-radius:0 0 10px 10px;box-shadow:0 8px 18px rgba(0,0,0,.04)}
        .result-block h5{color:#334155;border-bottom:1px solid var(--border)}
        .summary-block{background:#eef2f7;border:1px solid var(--border);border-radius:10px}
        .summary-block ul{list-style:none;margin:0;padding:0;display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:.6em 1em}
        .summary-block li{background:#fff;border:1px solid var(--border);padding:.5em .75em;border-radius:8px}
        /* Checklist visual polish and correct symbols */
        .checklist{margin:0}
        .checklist li{padding:.5em .75em;border-radius:6px;border:1px solid}
        /* Ensure no duplicate icons from pseudo-elements */
        .checklist li.pass::before, .checklist li.fail::before{content:'' !important; display:none}
        .checklist li.pass{background:#ecfdf3;border-color:#c6f6d5;color:#14532d}
        .checklist li.fail{background:#fef2f2;border-color:#fecaca;color:#7f1d1d}
        .checklist li .ic{display:inline-block;width:1.25em;text-align:center;font-weight:700;margin-right:.4em}
        /* Chips and controls */
        .chip{display:inline-block;padding:.15em .5em;border-radius:999px;font-size:.8em;border:1px solid var(--border);background:#fff;color:#334155;margin-left:.5em}
        .chip-ok{background:#ecfdf3;border-color:#bbf7d0;color:#166534}
        .chip-warn{background:#fffbeb;border-color:#fde68a;color:#92400e}
        .expand-controls{display:flex;gap:.5em;justify-content:flex-end;margin: .25em 0 1em}
        .expand-controls .secondary{background:#fff;color:#334155;border:1px solid var(--border);padding:.4em .8em;border-radius:6px;cursor:pointer}
        .expand-controls .secondary:hover{background:#f8fafc}
        /* Compact status chip for headers */
        .chip-status{display:inline-block;padding:.15em .55em;border-radius:999px;font-size:.82em;border:1px solid var(--border);font-weight:600;margin-left:.5em}
        .status-eligible{background:#14532d;color:#e6faef;border-color:#14532d}
        .status-likely-eligible{background:#0f4c5c;color:#e7f6fb;border-color:#0e7490}
        .status-unlikely-eligible{background:#fef3c7;color:#7c2d12;border-color:#f59e0b}
        .status-not-eligible{background:#7f1d1d;color:#fee2e2;border-color:#991b1b}
        .status-unable-to-assess{background:#334155;color:#e2e8f0;border-color:#334155}
        .status-error{background:#991b1b;color:#fee2e2;border-color:#991b1b}
        /* Professional-looking header colors */
        .eligible{background:linear-gradient(135deg,#15803d,#22c55e); color:#fff}
        .likely-eligible{background:linear-gradient(135deg,#0ea5e9,#06b6d4); color:#fff}
        .unlikely-eligible{background:linear-gradient(135deg,#f59e0b,#fbbf24); color:#111827}
        .not-eligible{background:linear-gradient(135deg,#dc2626,#ef4444); color:#fff}
        .unable-to-assess{background:linear-gradient(135deg,#64748b,#94a3b8); color:#fff}
        .error{background:linear-gradient(135deg,#b91c1c,#ef4444); color:#fff}
        /* Dark mode overrides */
        body.dark{--bg:#0b1220;--card:#0f172a;--text:#e5e7eb;--muted:#94a3b8;--border:#1f2937;--primary:#60a5fa;--primary-dark:#3b82f6}
        body.dark{color:var(--text)}
        body.dark a{color:#93c5fd}
        body.dark nav .links a, body.dark .theme-toggle{color:#e5e7eb}
        body.dark .container{box-shadow:0 6px 18px rgba(0,0,0,.5)}
        body.dark .summary-block{background:#111827;border-color:#1f2937}
        body.dark .summary-block li{background:#0b1324;border-color:#1f2937;color:var(--text)}
        body.dark .result-block{border-color:#1f2937;background:#0b1324;box-shadow:0 8px 18px rgba(0,0,0,.5)}
        /* Dark mode form controls */
        body.dark input, body.dark select, body.dark textarea{background:#0b1324;border:1px solid #1f2937;color:var(--text)}
        body.dark .expand-controls .secondary{background:#0b1324;color:#e5e7eb;border:1px solid #1f2937}
        body.dark .expand-controls .secondary:hover{background:#111827}
        /* Dark mode checklist tints */
        body.dark .checklist li.pass{background:#0f2e1c;border-color:#14532d;color:#bbf7d0}
        body.dark .checklist li.fail{background:#2a1414;border-color:#7f1d1d;color:#fecaca}
        /* Dark mode PubMed cards */
        body.dark .pubmed-card{background:#0b1324;border-color:#1f2937}
        body.dark .pubmed-title{color:#e5e7eb}
        body.dark .pubmed-article-link{background:#0b1324;border-color:#1f2937;color:#93c5fd}
        body.dark .pubmed-article-link:hover{background:#111827}
        body.dark .pubmed-footer a{background:#0b1324;border-color:#1f2937;color:#93c5fd}
        body.dark .pubmed-footer a:hover{background:#111827}
        /* Dark mode API docs code blocks */
        body.dark pre{background:#0b1324;border-color:#1f2937;color:#e5e7eb}
        body.dark code{color:#e5e7eb}
        body.dark pre code{color:#e5e7eb}
        /* Dark mode navbar blends into background */
        body.dark nav { background-color: #0b1220; }
        /* Dark mode IGV inversion trick */
        body.dark #igv-container{background:#0b1324; filter: invert(0.92) hue-rotate(180deg)}
        /* Loader spinner */
        #loader{position:relative}
        #loader:after{content:'';display:inline-block;width:1em;height:1em;margin-left:.5em;border:2px solid #cfd4dc;border-top-color:var(--primary);border-radius:50%;animation:spin 1s linear infinite;vertical-align:-2px}
        @keyframes spin{to{transform:rotate(360deg)}}
        /* Dark background override */
        body.dark{ background: linear-gradient(180deg,#0b1220 0%, #0a1020 100%); }
        @media (max-width: 640px){ #assessment-form{grid-template-columns:1fr}}
    </style>
</head>
<body>
    <nav>
        <div class="links">
            <a href="/">Tool</a>
            <a href="/about">About/Methods</a>
            <a href="/cite">How to Cite</a>
            <a href="/api_docs">API</a>
        </div>
        <button id="theme-toggle" class="theme-toggle" aria-pressed="false" title="Toggle dark mode">🌙</button>
    </nav>
    <script>
        (function(){
            const key = 'theme';
            const btn = document.getElementById('theme-toggle');
            function apply(t){
                const dark = (t === 'dark');
                document.body.classList.toggle('dark', dark);
                if (btn) {
                    btn.textContent = dark ? '☀️' : '🌙';
                    btn.setAttribute('aria-pressed', String(dark));
                }
            }
            let saved = localStorage.getItem(key);
            if (!saved) saved = window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light';
            apply(saved);
            btn && btn.addEventListener('click', () => {
                const next = document.body.classList.contains('dark') ? 'light' : 'dark';
                localStorage.setItem(key, next);
                apply(next);
            });
        })();
    </script>
    <main class="container">
        {% block content %}{% endblock %}
        <p class="disclaimer">This tool is for informational purposes only and is not for clinical use. Results require manual verification.</p>
    </main>
</body>
</html>
    """
    index_html = """
{% extends "base.html" %}
{% block content %}
<style>
    .splice-prompt-buttons { 
        margin-top: 1em; 
        display: flex; 
        gap: 10px; 
    }
    .splice-prompt-buttons button {
        flex: 1;
        color:white;
        font-weight:bold;
        cursor:pointer; 
        border:none; 
        width:auto;
        padding:.75em 1.5em; 
        border-radius: 4px; 
        font-size: 0.9em;
    }
    .splice-prompt-buttons button.yes-btn { background-color: #28a745; }
    .splice-prompt-buttons button.yes-btn:hover { background-color: #218838; }
    .splice-prompt-buttons button.no-btn { background-color: #dc3545; }
    .splice-prompt-buttons button.no-btn:hover { background-color: #c82333; }
    .warning { color: #dc3545; font-weight: 600; }
</style>

<h3>AVEC: Automated Variant Eligibility Calculator</h3>
<p>Enter a variant to assess its eligibility for ASO therapy.</p>
<p class="warning">Results depend on underlying data sources and external services and may be incomplete or unavailable. This tool does not replace clinical judgement or a physician!</p>
<form id="assessment-form"> <label for="query">Variant:</label> <input id="query" required placeholder="e.g., NM_015427.4:c.1054G>A"> <button type="submit">Assess</button> </form>
<div id="loader">Assessing...</div>
<div id="results"></div>

<hr style="margin: 2em 0;">

<h4 id="batch-toggle" style="cursor: pointer; user-select: none;">
    Batch Processing &#9662;
</h4>

<div id="batch-content" style="display: none;">
    <p>Upload a .csv, .txt, or .xlsx file with one variant per line in the first column.</p>
    <form id="batch-form">
        <label for="batch-file">Batch File:</label>
        <input type="file" id="batch-file" name="file" accept=".csv,.txt,.xlsx,.tsv" required>
        <button type="submit">Process Batch</button>
    </form>
    <div id="batch-loader" style="display:none; text-align: center; padding: 1em;">
        Processing file... This may take several minutes for large files.
    </div>
</div>

<script>
// Persist user selections across reassessments within the page lifecycle
window.userSpliceInput = null; // 'yes' | 'no' | null
window.userMoaInput = null;    // 'GoF' | 'LoF' | null
document.getElementById('assessment-form').addEventListener('submit', async function(e) {
    e.preventDefault();
    const resultsDiv = document.getElementById('results');
    const loader = document.getElementById('loader');
    resultsDiv.style.display = 'none';
    resultsDiv.innerHTML = '';
    loader.style.display = 'block';
    // New query: reset stored user inputs
    window.userSpliceInput = null;
    window.userMoaInput = null;
    
    // Initial assessment payload contains only the query
    const payload = { 
        query: document.getElementById('query').value 
        // splice_user_input is omitted, so backend defaults to DB check
    }; 
    
    try {
        const response = await fetch('/assess', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(payload) });
        const data = await response.json();
        displayResults(data);
    } catch (error) {
        console.error("Fetch Error:", error);
        displayResults({ classification: "Error", reason: "Could not connect to the server." });
    } finally {
        loader.style.display = 'none';
    }
});

document.getElementById('batch-form').addEventListener('submit', async function(e) {
    e.preventDefault();
    const batchLoader = document.getElementById('batch-loader');
    const fileInput = document.getElementById('batch-file');

    if (fileInput.files.length === 0) {
        alert("Please select a file to upload.");
        return;
    }

    batchLoader.style.display = 'block';
    const formData = new FormData();
    formData.append('file', fileInput.files[0]);

    try {
        const response = await fetch('/batch_assess', {
            method: 'POST',
            body: formData
        });

        if (response.ok) {
            const blob = await response.blob();
            const url = window.URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.style.display = 'none';
            a.href = url;
            const timestamp = new Date().toISOString().slice(0, 19).replace(/[-T:]/g, "");
            a.download = `avec_batch_results_${timestamp}.xlsx`;
            document.body.appendChild(a);
            a.click();
            window.URL.revokeObjectURL(url);
            a.remove();
        } else {
            const errorData = await response.json();
            alert(`Error processing file: ${errorData.error}`);
        }
    } catch (error) {
        console.error("Batch fetch error:", error);
        alert("A critical error occurred while communicating with the server.");
    } finally {
        batchLoader.style.display = 'none';
        fileInput.value = '';
    }
});

document.getElementById('batch-toggle').addEventListener('click', function() {
    const content = document.getElementById('batch-content');
    const toggle = this;

    if (content.style.display === 'none') {
        content.style.display = 'block';
        toggle.innerHTML = 'Batch Processing &#9652;'; // Up arrow
    } else {
        content.style.display = 'none';
        toggle.innerHTML = 'Batch Processing &#9662;'; // Down arrow
    }
});

function renderIGV(containerId, data) {
    const container = document.getElementById(containerId);
    if (!container || !data || !data.locus) {
        console.error("IGV container or data is missing, cannot render viewer.");
        return;
    }
    
    container.style.display = 'block';
    container.innerHTML = '';
    const options = { genome: "hg38", locus: data.locus };

    igv.createBrowser(container, options)
        .then(function (browser) {
            const tracksToLoad = [];
            if (data.variantTrack && data.variantTrack.features) {
                tracksToLoad.push({ type: "annotation", name: data.variantTrack.name, features: data.variantTrack.features, color: "red", displayMode: "EXPANDED", height: 35 });
            }
            if (data.domainTrack && data.domainTrack.features) {
                tracksToLoad.push({ type: "annotation", name: data.domainTrack.name, features: data.domainTrack.features, color: "#D46A6A", displayMode: "EXPANDED", height: 35 });
            }
            if (tracksToLoad.length > 0) browser.loadTrackList(tracksToLoad);
        });
}
// --- HELPER FUNCTION FOR SPLICE VALIDATION ---
async function reassessWithSpliceInput(spliceInput) {
    const resultsDiv = document.getElementById('results');
    const loader = document.getElementById('loader');
    resultsDiv.style.display = 'none';
    resultsDiv.innerHTML = '';
    loader.style.display = 'block';

    // Create the payload with the new 'splice_user_input' flag
    const payload = { 
        query: document.getElementById('query').value,
        splice_user_input: spliceInput, // 'yes' or 'no'
        // include previously chosen MoA if available to avoid losing state
        moa_user_input: window.userMoaInput
    };

    try {
        const response = await fetch('/assess', { 
            method: 'POST', 
            headers: { 'Content-Type': 'application/json' }, 
            body: JSON.stringify(payload) 
        });
        const data = await response.json();

        // Add a note to the re-assessed data so the user knows what happened
        if (data.assessments && data.assessments.Splice_Switching) {
            const originalReason = data.assessments.Splice_Switching.reason || "";
            let prefix = (spliceInput === 'yes') ? 
                "<strong>Re-assessed based on user-provided 'Yes'.</strong>" : 
                "<strong>Re-assessed based on user-provided 'No'.</strong>";
            data.assessments.Splice_Switching.reason = `${prefix} ${originalReason}`;
        }
        // persist the user's splice selection
        window.userSpliceInput = spliceInput;

        displayResults(data);
    } catch (error) {
        console.error("Fetch Error:", error);
        displayResults({ classification: "Error", reason: "Could not connect to the server." });
    } finally {
        loader.style.display = 'none';
    }
}

// --- NEW HELPER FUNCTION FOR MoA SELECTION ---
async function reassessWithMoaInput(moaChoice) {
    const resultsDiv = document.getElementById('results');
    const loader = document.getElementById('loader');
    resultsDiv.style.display = 'none';
    resultsDiv.innerHTML = '';
    loader.style.display = 'block';

    const payload = {
        query: document.getElementById('query').value,
        moa_user_input: moaChoice, // 'GoF' or 'LoF'
        // include previously chosen splice input if available to avoid losing state
        splice_user_input: window.userSpliceInput
    };

    try {
        const response = await fetch('/assess', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(payload)
        });
        const data = await response.json();

        // Add a note in the summary to reflect the user override
        if (data.summary) {
            const chosen = moaChoice;
            if (!data.summary.note) data.summary.note = '';
            data.summary.note = (data.summary.note ? data.summary.note + ' ' : '') +
                `user selection: ${chosen}`;
        }

        // persist the user's MoA selection
        window.userMoaInput = moaChoice;

        displayResults(data);
    } catch (error) {
        console.error('Fetch Error:', error);
        displayResults({ classification: 'Error', reason: 'Could not connect to the server.' });
    } finally {
        loader.style.display = 'none';
    }
}

// --- NEW EVENT LISTENER (for 'results' div) ---
document.getElementById('results').addEventListener('click', async function(e) {
    // Check if the "Yes" button was clicked
    if (e.target && e.target.id === 'splice-validation-yes') {
        e.preventDefault();
        reassessWithSpliceInput('yes');
    }
    
    // Check if the "No" button was clicked
    if (e.target && e.target.id === 'splice-validation-no') {
        e.preventDefault();
        reassessWithSpliceInput('no');
    }
    // MoA selection buttons
    if (e.target && e.target.id === 'moa-choice-gof') {
        e.preventDefault();
        reassessWithMoaInput('GoF');
    }
    if (e.target && e.target.id === 'moa-choice-lof') {
        e.preventDefault();
        reassessWithMoaInput('LoF');
    }
});


// --- MODIFIED displayResults FUNCTION (NOW CLEANED) ---
function displayResults(data) {
    const resultsDiv = document.getElementById('results');
    resultsDiv.style.display = 'block';
    resultsDiv.innerHTML = '';

    if (!data || data.classification === "Error" || data.classification === "Unable to Assess") {
        const classificationClass = (data.classification || "Error").toLowerCase().replace(/ /g, '-');
        resultsDiv.innerHTML = `<div class="result-header ${classificationClass}"><h4>${data.classification || "Error"}</h4></div><div class="result-block"><p>${data.reason || "An unknown error occurred."}</p></div>`;
        return;
    }

    let html = '';\n\n    // Precompute MoA-related values for use across sections\n    const assessmentsObj = data.assessments || {};\n    const summaryObj = data.summary || {};\n    let moaList = summaryObj.moa || [];\n    let note = summaryObj.note || '';\n    let resolvedMoa = summaryObj.resolved_moa || window.userMoaInput || null;\n    if (!resolvedMoa) { const m = note.match(/user selection:\\s*(GoF|LoF)/i); if (m) resolvedMoa = m[1]; }\n    const hasN1C = !!(assessmentsObj && (assessmentsObj.N1C_Registry_Check || assessmentsObj.N1C_Assessed_Variants));

    if (data.summary) {
        let geneHTML = data.summary.gene || 'N/A';
        if (data.summary.gene_url) {
            geneHTML = `<a href="${data.summary.gene_url}" target="_blank" rel="noopener noreferrer">${geneHTML}</a>`;
        }
        
        let haploHTML = (data.summary.haploinsufficiency && data.summary.haploinsufficiency.text) || 'N/A';
        if (data.summary.haploinsufficiency && data.summary.haploinsufficiency.url) {
            haploHTML = `<a href="${data.summary.haploinsufficiency.url}" target="_blank" rel="noopener noreferrer">${haploHTML}</a>`;
        }
        // Transcript link to Ensembl transcript page when available
        const txId = data.summary.transcript_id || null;
        const transcriptHTML = txId
            ? (String(txId).startsWith('ENST')
                ? `<a href="https://www.ensembl.org/Homo_sapiens/Transcript/Summary?t=${encodeURIComponent(txId)}" target="_blank" rel="noopener noreferrer">${txId}</a>`
                : `<a href="https://www.ensembl.org/Homo_sapiens/Search/Results?q=${encodeURIComponent(txId)};site=ensembl" target="_blank" rel="noopener noreferrer">${txId}</a>`)
            : 'N/A';
        // Mode of Inheritance links to the same page as Gene (external shows both)
        let moiHTML = 'N/A';
        const moiArr = data.summary.moi || [];
        if (moiArr.length > 0) {
            if (data.summary.gene_url) {
                moiHTML = moiArr.map(m => `<a href="${data.summary.gene_url}" target="_blank" rel="noopener noreferrer">${m}</a>`).join(', ');
            } else {
                moiHTML = moiArr.join(', ');
            }
        }
        // Orphanet link for the gene
        const orphaURL = data.summary.gene ? `https://www.orpha.net/en/disease/gene/${encodeURIComponent(data.summary.gene)}` : null;
        const orphaLinkHTML = orphaURL ? `<a href="${orphaURL}" target="_blank" rel="noopener noreferrer">${data.summary.gene}</a>` : 'N/A';

        html += `<div class="summary-block"><h4>Query Summary</h4><ul>
                        <li><strong>Gene:</strong> ${geneHTML}</li>
                        <li><strong>Transcript:</strong> ${transcriptHTML}</li>
                        <li><strong>Mode of Inheritance:</strong> ${moiHTML}</li>
                        <li><strong>Haploinsufficiency:</strong> ${haploHTML}</li>
                        <li><strong>Orphanet:</strong> ${orphaLinkHTML}</li>
                        <li><strong>Known Mechanism:</strong> ${(data.summary.moa && data.summary.moa.join(', ')) || 'N/A'}</li>
                        ${data.summary.note ? `<li><strong>Note:</strong> ${data.summary.note}</li>` : ''}
                        </ul></div>`;
    }

    html += '<div id="igv-container" style="display:none;"></div>';
    
    // Mechanism of Action block below visualization
    if (data.summary) { moaList = (data.summary.moa || []);
        note = data.summary.note || '';
        
        resolvedMoa = data.summary.resolved_moa || window.userMoaInput || resolvedMoa;
        if (!resolvedMoa) { const m = note.match(/user selection:\\s*(GoF|LoF)/i); if (m) resolvedMoa = m[1]; }
        if (!hasN1C) {
            if (resolvedMoa) {
                html += `<div class="result-block">
                            <h5>Mechanism of Action</h5>
                            <p>User selected: ${resolvedMoa}</p>
                        </div>`;
            } else if (moaList.length > 0) {
                html += `<div class="result-block">
                            <h5>Mechanism of Action</h5>
                            <p>${moaList.join(', ')}</p>
                        </div>`;
            }
        }

        // Ensure MoA is shown when multiple or none are known (even if N1C links exist)
        const showMoaBlock = (moaList.length === 0 || moaList.length > 1);
        if (showMoaBlock && hasN1C) {
            if (resolvedMoa) {
                const context = (moaList && moaList.length > 1) ? ` <span class=\"note\">(Known: ${moaList.join(', ')})</span>` : '';
                html += `<div class=\"result-block\">\n                            <h5>Mechanism of Action</h5>\n                            <p>User selected: ${resolvedMoa}${context}</p>\n                        </div>`;
            } else {
                const body = (moaList && moaList.length > 1)
                    ? `Multiple mechanisms reported: ${moaList.join(', ')}`
                    : `No established mechanism of action is known.`;
                html += `<div class=\"result-block\">\n                            <h5>Mechanism of Action</h5>\n                            <p>${body} Please select the appropriate mechanism below.</p>\n                        </div>`;
            }
        }
        // Also show when no MoA is known and there is no N1C block
        if (moaList.length === 0 && !resolvedMoa && !hasN1C) {
            html += `<div class=\"result-block\">\n                        <h5>Mechanism of Action</h5>\n                        <p>No established mechanism of action is known. Please select the appropriate mechanism below.</p>\n                    </div>`;
        }
    }
    
    // Literature accordion (PubMed and Google Scholar) above assessments
    if (data.summary && data.summary.gene) {
        const gene = data.summary.gene;
        const pubmedQuery = `(${gene}) AND (ASO OR AON OR \"Antisense oligonucleotide\")`;
        const pubmedURL = 'https://pubmed.ncbi.nlm.nih.gov/?term=' + encodeURIComponent(pubmedQuery);
        // const scholarQuery = `(${gene}) AND (ASO OR AON OR \"Antisense oligonucleotide\")`;
        // const scholarURL = 'https://scholar.google.com/scholar?q=' + encodeURIComponent(scholarQuery);

        html += `<div class="strategy-block">
                <div id="toggle-literature" class="result-header unable-to-assess accordion-toggle" style="cursor: pointer; user-select: none;">
                    <h4>Literature: PubMed <span class="chevron">&#9662;</span></h4>
                </div>
                <div id="content-literature" class="result-block" style="display: none;">
                    <p><strong>Query:</strong> ${pubmedQuery}</p>
                    <div id="pubmed-results" data-gene="${gene}">Loading PubMed results...</div>
                    <div class="pubmed-footer"><a href="${pubmedURL}" target="_blank" rel="noopener noreferrer">Open full results on PubMed</a></div>
                    <p class="note" style="margin-top:0.5em;">Powered by NCBI E-utilities.</p>
                </div>
            </div>`;
    }

    if (data.assessments) {
        html += '<h4>Therapeutic Assessments</h4>';
        html += '<div class="expand-controls"><button id="expand-all" class="secondary">Expand All</button><button id="collapse-all" class="secondary">Collapse All</button></div>';
        // MoA selection block inside Therapeutic Assessments when needed
        if (!resolvedMoa && (moaList.length === 0 || moaList.length > 1)) {
            html += `<div class="strategy-block">`;
            html += `<div id="toggle-moa" class="result-header unable-to-assess accordion-toggle" style="cursor: pointer; user-select: none;">
                         <h4>Mechanism of Action: Selection Required <span class="chevron">&#9662;</span></h4>
                     </div>`;
            html += `<div id="content-moa" class="result-block" style="display: block;">
                        <p style="font-weight: bold;">There is no single known mechanism of action available. Please select the mechanism relevant to your variant:</p>
                        <div class="splice-prompt-buttons">
                            <button id="moa-choice-gof" type="button" class="yes-btn">GoF</button>
                            <button id="moa-choice-lof" type="button" class="no-btn">LoF</button>
                        </div>
                    `;
            // Insert contextual PubMed search for MoA right below the buttons
            (function(){
                const gene = (data.summary && data.summary.gene) ? data.summary.gene : '';
                const origQ = window.lastQuery || (document.getElementById('query') ? document.getElementById('query').value : '');
                const cdnaMatch = origQ ? origQ.match(/c\.[^\s\)\]]+/i) : null;
                const cdna = cdnaMatch ? cdnaMatch[0] : '';
                const refseqMatch = origQ ? origQ.match(/(NM_[0-9]+\.?[0-9]*)/i) : null;
                const refseq = refseqMatch ? refseqMatch[1] : '';
                const parts = [];
                if (gene && cdna) parts.push(`(${gene} AND ${cdna})`);
                if (refseq) parts.push(refseq);
                const prefix = parts.length ? `(${parts.join(' OR ')})` : '';
                const terms = '(gain-of-function OR GoF OR loss-of-function OR LoF OR "Gain of function" OR "Loss of function" OR "Mechanism of action")';
                const moaQuery = prefix ? `${prefix} AND ${terms}` : terms;
                const pmUrl = 'https://pubmed.ncbi.nlm.nih.gov/?term=' + encodeURIComponent(moaQuery);
                html += `\n                    <div class=\"pubmed-supplement\">\n                        <p><strong>Related PubMed Search:</strong> ${moaQuery}</p>\n                        <div id=\"pubmed-moa-results\" data-query=\"${moaQuery}\">Loading PubMed results...</div>\n                        <div class=\"pubmed-footer\"><a href=\"${pmUrl}\" target=\"_blank\" rel=\"noopener noreferrer\">Open full results on PubMed</a></div>\n                        <p class=\"note\" style=\"margin-top:0.5em;\">Powered by NCBI E-utilities.</p>\n                    </div>`;
            })();
            html += `</div></div>`;
        }
        for (const [strategy, result] of Object.entries(data.assessments)) {
            const strategyName = strategy.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase());

            // --- START: MODIFIED LOGIC FOR SPLICE_SWITCHING ---
            if (strategy === 'Splice_Switching' && result.user_validation_prompt === true) {
                // Render the special "Awaiting Validation" block with "Yes/No" buttons
                html += `<div class="strategy-block">`;
                html += `<div id="toggle-${strategy}" class="result-header unable-to-assess accordion-toggle" style="cursor: pointer; user-select: none;">
                             <h4>${strategyName}: Awaiting Validation <span class="chevron">&#9662;</span></h4>
                         </div>`;
                // Set to 'block' so the user sees the question immediately
                html += `<div id="content-${strategy}" class="result-block" style="display: block;">`; 
                // The new question
                html += `<p style="font-weight: bold;">Splicing effect was not found in databases. Is there a known splice-altering effect validated with qPCR or Transcriptomics?</p>`;
                // The new buttons
                html += `<div class="splice-prompt-buttons">
                             <button id="splice-validation-yes" type="button" class="yes-btn">Yes</button>
                             <button id="splice-validation-no" type="button" class="no-btn">No</button>
                         </div>`;
                // Insert contextual PubMed below buttons, then close block
                (function(){
                    const gene = (data.summary && data.summary.gene) ? data.summary.gene : '';
                    const origQ = window.lastQuery || (document.getElementById('query') ? document.getElementById('query').value : '');
                    const cdnaMatch = origQ ? origQ.match(/c\.[^\s\)\]]+/i) : null;
                    const cdna = cdnaMatch ? cdnaMatch[0] : '';
                    const refseqMatch = origQ ? origQ.match(/(NM_[0-9]+\.?[0-9]*)/i) : null;
                    const refseq = refseqMatch ? refseqMatch[1] : '';
                    const parts = [];
                    if (gene && cdna) parts.push(`(${gene} AND ${cdna})`);
                    if (refseq) parts.push(refseq);
                    const prefix = parts.length ? `(${parts.join(' OR ')})` : '';
                    const terms = '(Splicing OR Splice-altering OR \"intron retention\" OR \"splice-switch\" OR \"cryptic splice site activation\")';
                    const spliceQuery = prefix ? `${prefix} AND ${terms}` : terms;
                    const pmUrl = 'https://pubmed.ncbi.nlm.nih.gov/?term=' + encodeURIComponent(spliceQuery);
                    html += `\n                    <div class=\"pubmed-supplement\">\n                        <p><strong>Related PubMed Search:</strong> ${spliceQuery}</p>\n                        <div id=\"pubmed-splice-results\" data-query=\"${spliceQuery}\">Loading PubMed results...</div>\n                        <div class=\"pubmed-footer\"><a href=\"${pmUrl}\" target=\"_blank\" rel=\"noopener noreferrer\">Open full results on PubMed</a></div>\n                        <p class=\"note\" style=\"margin-top:0.5em;\">Powered by NCBI E-utilities.</p>\n                    </div>`;
                })();
                html += `</div></div>`;
            
            } else {
                // This is the "else" block: it contains all the *original* rendering logic
                const classificationClass = (result.classification || "Unable to Assess").toLowerCase().replace(/ /g, '-');
                
                html += `<div class="strategy-block">`;
                html += `<div id="toggle-${strategy}" class="result-header ${classificationClass} accordion-toggle" style="cursor: pointer; user-select: none;">
                             <h4>${strategyName} <span class="chip-status status-${classificationClass}">${result.classification || "N/A"}</span> <span class="chevron">&#9662;</span></h4>
                         </div>`;
                // Set to 'none' by default
                html += `<div id="content-${strategy}" class="result-block" style="display: none;">`;
                
                if (strategy === 'Exon_Skipping' && result.total_exon_number && result.gene_id && result.transcript_id) {
                    const ensemblLink = `https://www.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=${result.gene_id};t=${result.transcript_id}`;
                    html += `<p><strong>Target:</strong> <a href="${ensemblLink}" target="_blank" rel="noopener noreferrer">Exon ${result.total_exon_number}</a></p>`;
                }
                if(result.reason) html += `<p><strong>Reason:</strong> ${result.reason}</p>`;
                
                if(result.details) {
                    html += '<h5>Details</h5><ul>';
                    for (const [key, value] of Object.entries(result.details)) {
                        if (value && (value.startsWith('http') || value.startsWith('https'))) {
                            html += `<li><strong>${key}:</strong> <a href="${value}" target="_blank" rel="noopener noreferrer">${value}</a></li>`;
                        } else {
                            html += `<li><strong>${key}:</strong> ${value}</li>`;
                        }
                    }
                    html += '</ul>';
                }

                if (result.checks && strategy !== 'Allele_Specific_Knockdown') {
                    const totalChecks = Object.keys(result.checks).length;
                    const passedChecks = Object.values(result.checks).filter(Boolean).length;
                    html += `<h5>Guideline Checks <span class="chip ${passedChecks===totalChecks ? 'chip-ok':'chip-warn'}">${passedChecks}/${totalChecks} passed</span></h5><ul class="checklist">`;
                    for (const [check, passed] of Object.entries(result.checks)) {
                        html += `<li class="${passed ? 'pass' : 'fail'}">${check}</li>`;
                    }
                    html += '</ul>';
                }
                
                if (result.pathogenic_variant_counts) {
                    let evidenceHeader = '<h5>Evidence from Databases</h5>';
                    if (result.clinvar_url) {
                        evidenceHeader = `<h5><a href="${result.clinvar_url}" target="_blank" rel="noopener noreferrer">Evidence from Databases (ClinVar)</a></h5>`;
                    }
                    html += evidenceHeader;

                    let domainHTML = result.domain_count || 'N/A';
                    if (result.domain_names && result.domain_names.length > 0) {
                        domainHTML = result.domain_names.join(', ');
                    }

                    html += `<ul>
                                <li><strong>Fraction of Protein:</strong> ${result.frac_cds || 'N/A'}</li>
                                <li><strong>Overlapping Protein Domains:</strong> ${domainHTML}</li>
                                <li><strong>Pathogenic Variants in Exon:</strong><ul>
                                    <li>Missense: ${result.pathogenic_variant_counts.missense}</li>
                                    <li>Nonsense: ${result.pathogenic_variant_counts.nonsense}</li>
                                    <li>Frameshift: ${result.pathogenic_variant_counts.frameshift}</li>
                                    <li>In-frame Deletions: ${result.pathogenic_variant_counts.inframe_del}</li>
                                    <li>Splice Site: ${result.pathogenic_variant_counts.splice}</li>
                                </ul></li>
                            </ul>`;
                }
                html += `</div></div>`;
            }
            // --- END: MODIFIED LOGIC ---
        }
    }
    
    resultsDiv.innerHTML = html;\n\n    // Expand/Collapse all wiring\n    (function(){\n        const expandBtn = document.getElementById('expand-all');\n        const collapseBtn = document.getElementById('collapse-all');\n        if (expandBtn) expandBtn.addEventListener('click', () => {\n            document.querySelectorAll('.result-block').forEach(el => el.style.display='block');\n            document.querySelectorAll('.accordion-toggle').forEach(t => t.classList.add('open'));\n        });\n        if (collapseBtn) collapseBtn.addEventListener('click', () => {\n            document.querySelectorAll('.result-block').forEach(el => el.style.display='none');\n            document.querySelectorAll('.accordion-toggle').forEach(t => t.classList.remove('open'));\n        });\n    })();

    // Fetch PubMed articles via NCBI E-utilities and render inline
    (async () => {
        const container = document.getElementById('pubmed-results');
        if (!container) return;
        const gene = container.getAttribute('data-gene');
        const query = `(${gene}) AND (ASO OR AON OR "Antisense oligonucleotide")`;
        try {
            const esearchUrl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&retmax=20&term=' + encodeURIComponent(query);
            const esResp = await fetch(esearchUrl);
            const esJson = await esResp.json();
            const ids = (esJson.esearchresult && esJson.esearchresult.idlist) || [];
            if (!ids.length) { container.innerHTML = 'No PubMed results found.'; return; }
            const esummaryUrl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=' + ids.join(',');
            const sumResp = await fetch(esummaryUrl);
            const sumJson = await sumResp.json();
            const result = sumJson.result || {};
            const uids = result.uids || ids;
            let out = '<ul class="pubmed-list">';
            for (const uid of uids) {
                const rec = result[uid];
                if (!rec) continue;
                const title = rec.title || '(no title)';
                const journal = rec.fulljournalname || rec.source || '';
                const pubdate = rec.pubdate || rec.epubdate || '';
                const authors = Array.isArray(rec.authors) ? rec.authors.slice(0,3).map(a => a.name).join(', ') : '';
                const link = `https://pubmed.ncbi.nlm.nih.gov/${uid}/`;
                out += `<li class="pubmed-card">` +
                       `<div class="pubmed-title">${title}</div>` +
                       `<div class="pubmed-meta">${journal}${journal && pubdate ? ' � ' : ''}${pubdate}${authors ? ' � ' + authors : ''}</div>` +
                       `<div class="pubmed-card-footer"><a class="pubmed-article-link" href="${link}" target="_blank" rel="noopener noreferrer">View on PubMed</a></div>` +
                       `</li>`;
            }
            out += '</ul>';
            container.innerHTML = out;
        } catch (e) {
            container.innerHTML = 'Unable to load PubMed results. Use the link above.';
        }
    })();

    // Render contextual PubMed results for MoA/Splice queries if present
    (async () => {
        const targets = document.querySelectorAll('#pubmed-moa-results, #pubmed-splice-results');
        if (!targets || targets.length === 0) return;
        for (const container of targets) {
            const q = container.getAttribute('data-query');
            if (!q) continue;
            try {
                const esearchUrl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&retmax=20&term=' + encodeURIComponent(q);
                const esResp = await fetch(esearchUrl);
                const esJson = await esResp.json();
                const ids = (esJson.esearchresult && esJson.esearchresult.idlist) || [];
                if (!ids.length) { container.innerHTML = 'No PubMed results found.'; continue; }
                const esummaryUrl = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&retmode=json&id=' + ids.join(',');
                const sumResp = await fetch(esummaryUrl);
                const sumJson = await sumResp.json();
                const result = sumJson.result || {};
                const uids = result.uids || ids;
                let out = '<ul class=\"pubmed-list\">';
                for (const uid of uids) {
                    const rec = result[uid];
                    if (!rec) continue;
                    const title = rec.title || '(no title)';
                    const journal = rec.fulljournalname || rec.source || '';
                    const pubdate = rec.pubdate || rec.epubdate || '';
                    const authors = Array.isArray(rec.authors) ? rec.authors.slice(0,3).map(a => a.name).join(', ') : '';
                    const link = `https://pubmed.ncbi.nlm.nih.gov/${uid}/`;
                    out += `<li class=\"pubmed-card\">` +
                           `<div class=\"pubmed-title\">${title}</div>` +
                           `<div class=\"pubmed-meta\">${journal}${journal && pubdate ? ' ? ' : ''}${pubdate}${authors ? ' ? ' + authors : ''}</div>` +
                           `<div class=\"pubmed-card-footer\"><a class=\"pubmed-article-link\" href=\"${link}\" target=\"_blank\" rel=\"noopener noreferrer\">View on PubMed</a></div>` +
                           `</li>`;
                }
                out += '</ul>';
                container.innerHTML = out;
            } catch (e) {
                container.innerHTML = 'Unable to load PubMed results. Use the link above.';
            }
        }
    })();

    document.querySelectorAll('.accordion-toggle').forEach(toggle => {
        toggle.addEventListener('click', function() {
            const contentId = this.id.replace('toggle-', 'content-');
            const content = document.getElementById(contentId);
            const isHidden = content.style.display === 'none'; content.style.display = isHidden ? 'block' : 'none'; this.classList.toggle('open', isHidden);
        });
    });

    if (data.visualization) {
        renderIGV('igv-container', data.visualization);
    }
}
</script>
{% endblock %}
"""

    about_html = """
{% extends "base.html" %}
{% block content %}
<h3>About & Methods</h3>
<h4>Purpose</h4>
<p>AVEC (Automated Variant Eligibility Calculator) is a research tool designed to assess the eligibility of a variant for therapeutic antisense oligonucleotides (ASOs). It automates the analysis of key criteria to predict whether ASO therapy (splice-correcting, WT-upregulation, knockdown or exon skipping) is likely to restore a functional protein product, and assumes variants are nonsense, missense small indels or frameshift variants. For other variants please refer to the <a href="https://shorturl.at/YqphL" target="_blank" rel="noopener noreferrer"><b>guidelines</b></a>. If the underlying datasources are false or incomplete it is possible that the results are also false. Therefore we strongly advise manual validation of the results</p>

<h4>Methodology</h4>
<p>This is a prototype</p>
{% endblock %}
"""

    cite_html = """
{% extends "base.html" %}
{% block content %}
<h3>prototype
{% endblock %}
"""
    api_docs_html = """
{% extends "base.html" %}
{% block content %}
<h3>AVEC API Documentation</h3>
<p>
    The AVEC API provides programmatic access to the variant assessment tool. 
    You can retrieve the same detailed analysis available through the batch processing tool via a simple GET request.
</p>

<h4>Endpoint</h4>
<p>The base URL for the API endpoint is:</p>
<pre><code>{{ url_for('api_assess', _external=True) }}</code></pre>

<h4>Request</h4>
<p>The API accepts <strong>GET</strong> requests with a single required query parameter.</p>
<ul>
    <li><strong>Parameter:</strong> <code>query</code></li>
    <li><strong>Description:</strong> The variant to assess in a recognized HGVS-like format.</li>
    <li><strong>Examples:</strong> <code>NM_015427.4:c.1054G>A</code>, <code>FKTN c.1312G>A</code></li>
</ul>

<h4>Example Usage (cURL)</h4>
<pre><code>curl -X GET "{{ url_for('api_assess', _external=True) }}?query=NM_000552.4:c.545G>A"</code></pre>

<h4>Response</h4>
<p>The API returns a JSON object containing the full assessment, structured identically to the data used by the web interface.</p>
<ul>
    <li>On success (HTTP 200), the response will contain `summary` and `assessments` objects.</li>
    <li>On failure (e.g., invalid query or server error), it will return a JSON object with an `error` key.</li>
</ul>

<h5>Example Successful Response Snippet</h5>
<pre><code>{
  "assessments": {
    "Exon_Skipping": {
      "checks": {
        "Benign splice variant found": false,
        "Is <10% of Protein": true,
        "Is In-Frame": true,
        ...
      },
      "classification": "Likely Eligible",
      "domain_count": 0,
      "frac_cds": "3.55%",
      "reason": "Exon meets the primary criteria for a skippable exon."
    }
  },
  "summary": {
    "gene": "DMD",
    "haploinsufficiency": {
      "text": "No evidence",
      "url": "https://search.clinicalgenome.org/kb/genes/HGNC:2928"
    },
    "moa": [],
    "moi": [ "X-linked" ],
    "transcript_id": "ENST00000357033.9"
  },
  "visualization": { ... }
}</code></pre>

<h4>Fair Use</h4>
<p>This is a free and open research tool. Please limit requests to a reasonable rate to ensure service availability for all users. For very large batch jobs, we recommend using the file upload feature on the main page.</p>

{% endblock %}
"""

    # Write files with explicit UTF-8 encoding
    with open('templates/base.html', 'w', encoding='utf-8') as f: f.write(base_html)
    with open('templates/index.html', 'w', encoding='utf-8') as f: f.write(index_html)
    with open('templates/about.html', 'w', encoding='utf-8') as f: f.write(about_html)
    with open('templates/cite.html', 'w', encoding='utf-8') as f: f.write(cite_html)
    with open('templates/api_docs.html', 'w', encoding='utf-8') as f: f.write(api_docs_html)
# Run the setup function to ensure templates exist
setup_templates()

# ===== CONFIGURATION =====
ENSEMBL_REST = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json", "Accept": "application/json"}
N1C_API_URL = "https://gene-registry.onrender.com/api/data?table=N1C_projects" 
N1C_API_ASSESSED_URL = "https://gene-registry.onrender.com/api/data?table=assessed_variants"
# Global DataFrames to be loaded at startup
clingen_df: Optional[pd.DataFrame] = None
goflof_df: Optional[pd.DataFrame] = None
splicevar_df: Optional[pd.DataFrame] = None
n1c_variants_df: Optional[pd.DataFrame] = None 
n1c_assessed_df: Optional[pd.DataFrame] = None 
sscvdb_df: Optional[pd.DataFrame] = None 

# Supplementary N1C table with gene-level features (uORF, NAT, PE)
n1c_supp_df: Optional[pd.DataFrame] = None 

# --- Data Loading ---
def load_databases():
    """
    Loads all necessary data files and fetches N1C registry data
    into global pandas DataFrames.
    """
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    DATA_DIR = os.path.join(BASE_DIR, 'data')

    clingen_path = os.path.join(DATA_DIR, 'Clingen-Curation-Activity-Summary-Report-2025-10-15.csv')
    goflof_path = os.path.join(DATA_DIR, 'goflof_HGMD2019_v032021_allfeat.csv')
    splicevar_path = os.path.join(DATA_DIR, 'splicevardb.xlsx')
    sscvdb_path = os.path.join(DATA_DIR, 'SSCVDB.xlsx')
    n1c_supp_path = os.path.join(DATA_DIR, 'N1C_Variant_Supp_Table.xlsx')

    global clingen_df, goflof_df, splicevar_df, n1c_variants_df, n1c_assessed_df, sscvdb_df, n1c_supp_df
    print("Loading databases...")
    try:
        clingen_df = pd.read_csv(clingen_path).set_index('gene_symbol')
        goflof_df = pd.read_csv(goflof_path).set_index('GENE')
        
        # Load SpliceVarDB from Excel
        splicevar_df = pd.read_excel(splicevar_path)
        # Sanitize SpliceVarDB data (critical for lookups)
        splicevar_df.columns = splicevar_df.columns.str.strip()
        if 'hgvs' in splicevar_df.columns and 'gene' in splicevar_df.columns:
            splicevar_df['hgvs'] = splicevar_df['hgvs'].astype(str).str.strip()
            splicevar_df['gene'] = splicevar_df['gene'].astype(str).str.strip()
        print(f"Loaded and sanitized SpliceVarDB data for {len(splicevar_df)} variants.")

        # Load SSCVDB from Excel
        try:
            sscvdb_df = pd.read_excel(sscvdb_path)
            sscvdb_df.columns = sscvdb_df.columns.str.strip()
            if 'Variant ID' in sscvdb_df.columns:
                sscvdb_df['Variant ID'] = sscvdb_df['Variant ID'].astype(str).str.strip()
            print(f"Loaded SSCVDB entries: {len(sscvdb_df)}")
        except Exception as e:
            print(f"Warning: Could not load SSCVDB.xlsx: {e}")

        # Load N1C supplementary table (uORF / NAT / PE per gene)
        try:
            n1c_supp_df = pd.read_excel(n1c_supp_path)
            n1c_supp_df.columns = n1c_supp_df.columns.str.strip()
            for col in ['Gene', 'uORF', 'NAT', 'PE']:
                if col in n1c_supp_df.columns:
                    n1c_supp_df[col] = n1c_supp_df[col].astype(str).str.strip()
            print(f"Loaded N1C supplementary table entries: {len(n1c_supp_df)}")
        except Exception as e:
            print(f"Warning: Could not load N1C_Variant_Supp_Table.xlsx: {e}")

        # Fetch and load N1C N-of-1 Projects data
        print("Fetching N=1 Collaborative Projects registry data...")
        response = requests.get(N1C_API_URL, timeout=30)
        response.raise_for_status() # Will raise an error if the request fails
        n1c_data = response.json()
        n1c_variants_df = pd.DataFrame(n1c_data)
        
        # Sanitize the N1C columns we will search on
        if 'Gene' in n1c_variants_df.columns:
            n1c_variants_df['Gene'] = n1c_variants_df['Gene'].astype(str).str.strip()
        if 'Coding DNA change (c.)' in n1c_variants_df.columns:
            n1c_variants_df['Coding DNA change (c.)'] = n1c_variants_df['Coding DNA change (c.)'].astype(str).str.strip()
        print(f"Successfully loaded {len(n1c_variants_df)} projects from the N1C registry.")

        # Fetch and load N1C assessed variants (curated) data
        print("Fetching N1C assessed variants data...")
        response2 = requests.get(N1C_API_ASSESSED_URL, timeout=30)
        response2.raise_for_status()
        n1c_assessed_data = response2.json()
        n1c_assessed_df = pd.DataFrame(n1c_assessed_data)
        # Sanitize columns commonly used for matching
        if 'Gene' in n1c_assessed_df.columns:
            n1c_assessed_df['Gene'] = n1c_assessed_df['Gene'].astype(str).str.strip()
        # Normalize any potential c. notation columns to a unified helper accessor later
        for col in list(n1c_assessed_df.columns):
            if isinstance(col, str):
                n1c_assessed_df[col] = n1c_assessed_df[col].astype(str).str.strip()
        print(f"Successfully loaded {len(n1c_assessed_df)} curated assessed variants from the N1C registry.")

    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"An error occurred during database loading: {e}")
        exit(1)

def _format_sscvdb_variant_id_from_vep(vep_entry: Dict[str, Any]) -> Optional[str]:
    """Formats a VEP entry to SSCVDB Variant ID style: chr<chrom>-<pos>-<ref>-<alt>.
    Returns None if required fields are missing or allele string is ambiguous."""
    try:
        chrom = str(vep_entry.get('seq_region_name', '')).strip()
        pos = vep_entry.get('start')
        allele = vep_entry.get('allele_string')
        if not chrom or not pos or not allele:
            return None
        if '/' not in allele:
            return None
        ref, alt = allele.split('/', 1)
        chrom_prefixed = f"chr{chrom}" if not chrom.lower().startswith('chr') else chrom
        return f"{chrom_prefixed}-{pos}-{ref}-{alt}"
    except Exception:
        return None

class EnsemblClient:
    def __init__(self, base_url=ENSEMBL_REST, headers=HEADERS, delay=0.1):
        self.base_url = base_url.rstrip('/')
        self.session = requests.Session()
        self.session.headers.update(headers)
        self.delay = delay

    def _get(self, path, params=None, max_retries=5):
        url = f"{self.base_url}{path}"
        backoff = 1.0
        for attempt in range(max_retries):
            time.sleep(self.delay)
            try:
                resp = self.session.get(url, params=params, timeout=30)
                if resp.status_code == 200:
                    try: return resp.json()
                    except ValueError: return resp.text
                elif resp.status_code in (429, 503):
                    wait = float(resp.headers.get('Retry-After', backoff))
                    time.sleep(wait); backoff *= 2
                elif 500 <= resp.status_code < 600:
                    time.sleep(backoff); backoff *= 2
                else:
                    if 400 <= resp.status_code < 500: return None
            except requests.RequestException:
                if attempt + 1 == max_retries: raise
                time.sleep(backoff); backoff *= 2
        return None

    def lookup_id_expand(self, identifier): return self._get(f"/lookup/id/{identifier}", params={'expand': '1'})
    def vep_hgvs(self, hgvs_string): return self._get(f"/vep/human/hgvs/{hgvs_string.strip()}", params={'variant_class': 1})
    def get_cds_sequence(self, transcript_id):
        data = self._get(f"/sequence/id/{transcript_id}", params={"type": "cds"})
        return data.get("seq") if isinstance(data, dict) else None
    def get_domains(self, protein_id):
        all_features = self._get(f"/overlap/translation/{protein_id}", params={"feature": "protein_feature"})
        if not all_features or not isinstance(all_features, list): return []
        domain_sources = {'CDD','Pfam','SMART','PROSITE profiles','PROSITE patterns','SUPERFAMILY','PRINTS','TIGRFAM','ProDom'}
        preliminary_domains = [f for f in all_features if f.get('type') in domain_sources]
        unique_interpro_domains = {f['interpro']: f for f in preliminary_domains if f.get('interpro')}
        return list(unique_interpro_domains.values())
    def overlap_region_variation(self, chrom, start, end):
        data = self._get(f"/overlap/region/human/{chrom}:{start}-{end}", params={'feature': 'variation'})
        return data if isinstance(data, list) else []
    def get_overlapping_genes(self, gene_id):
        """Fetches all genes that overlap with a given Ensembl Gene ID."""
        data = self._get(f"/overlap/id/{gene_id}", params={"feature": "gene"})
        return data if isinstance(data, list) else []
    def lookup_symbol(self, symbol):
        """Fetches gene data for a given symbol."""
        data = self._get(f"/lookup/symbol/human/{symbol}", params={'expand': '0'})
        return data if isinstance(data, dict) else None
    
# --- Helper & Parsing Functions ---

def _evaluate_splice_variant_position(variant_hgvs: str, vep_data: Dict[str, Any], details: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    """
    Assesses a validated splice-altering variant based on its genomic position.
    This logic is shared by database-found variants and user-validated variants.
    """
    core_hgvs_match = re.search(r'(c\..*)', variant_hgvs, re.IGNORECASE)
    if not core_hgvs_match: 
        return None
    core_canonical_hgvs = core_hgvs_match.group(1).lower()

    result = {"details": details}
    consequence_terms = set(vep_data.get('transcript_consequences', [{}])[0].get('consequence_terms', []))
    
    is_intronic_by_consequence = 'intron_variant' in consequence_terms or 'splice_acceptor_variant' in consequence_terms or 'splice_donor_variant' in consequence_terms
    is_intronic_by_notation = '+' in core_canonical_hgvs or '-' in core_canonical_hgvs
    is_intronic = is_intronic_by_consequence or is_intronic_by_notation

    if is_intronic:
        dist_match = re.search(r'[+-](\d+)', core_canonical_hgvs)
        if not dist_match: 
            return None 
        dist = int(dist_match.group(1))
        
        if '+' in core_canonical_hgvs: # Downstream (e.g., c.123+1G>A)
            if dist <= 5:
                result.update({"classification": "Not Eligible", "reason": "Variant is a validated splice-altering variant located too close to the canonical splice site (<=+5bp)."})
            elif 6 <= dist <= 50:
                result.update({"classification": "Unlikely Eligible", "reason": "Variant is a validated splice-altering variant located near the canonical splice site (+6-+50bp)."})
            else:
                result.update({"classification": "Likely Eligible", "reason": "Variant is a validated splice-altering variant in a favorable deep-intronic position (>+50bp)."})
        
        elif '-' in core_canonical_hgvs: # Upstream (e.g., c.124-2A>G)
            if dist <= 5:
                result.update({"classification": "Not Eligible", "reason": "Variant is a validated splice-altering variant located too close to the canonical splice site (>=-5bp)."})
            elif 6 <= dist <= 100:
                result.update({"classification": "Unlikely Eligible", "reason": "Variant is a validated splice-altering variant located near the canonical splice site (-6-(-100b)p)."})
            else:
                result.update({"classification": "Likely Eligible", "reason": "Variant is a validated splice-altering variant in a favorable deep-intronic position (<-100bp)."})
        
    elif any(c in consequence_terms for c in ['missense_variant', 'synonymous_variant']):
        if 'splice_region_variant' in consequence_terms:
            result.update({"classification": "Not Eligible", "reason": "This validated exonic splice-altering variant is within the canonical splice region, making it high-risk."})
        else:
            result.update({"classification": "Likely Eligible", "reason": "This validated exonic splice-altering variant is outside the immediate splice region, making it a potential candidate for correction."})
    
    else:
        # Fallback if it's splice-altering but not in a recognized position
        result.update({"classification": "Unlikely Eligible", "reason": "Variant is a validated splice-altering variant, but also presumed to cause other effects (e.g. it is a nonsense variant)."})

    return result

def parse_hgvs_query(query: str) -> Tuple[Optional[str], Optional[str]]:
    """
    Parses a query string into a VEP-compatible HGVS string and an optional gene symbol.
    """
    query = query.strip()
    
    # Pattern 1: Handles formats like "NM_12345.6:c.123A>G" or "GENE:c.123A>G"
    match_colon = re.search(r'([^:]+):([cgnmp]\..*)', query, re.IGNORECASE)
    if match_colon:
        identifier = match_colon.group(1).strip()
        variant = match_colon.group(2).strip()
        hgvs_string = f"{identifier}:{variant}"
        # If the identifier is a transcript, we don't have a gene symbol from the query
        if identifier.startswith("NM_") or identifier.startswith("ENST"):
            return hgvs_string, None
        else:
            return hgvs_string, identifier

    # Pattern 2: Handles formats like "GENE c.123A>G"
    match_space = re.search(r'([A-Z0-9\-_]+)\s+([cgnmp]\..*)', query, re.IGNORECASE)
    if match_space:
        gene = match_space.group(1).strip()
        variant = match_space.group(2).strip()
        return f"{gene}:{variant}", gene
        
    return None, None

def classify_variant_clinsig(clinsig_field):
    if clinsig_field is None: return 'other'
    vals = [v.lower() for v in (clinsig_field if isinstance(clinsig_field, list) else [clinsig_field]) if isinstance(v, str)]
    if any('pathogenic' in v for v in vals) and not any('likely' in v for v in vals): return 'pathogenic'
    if any('likely pathogenic' in v for v in vals): return 'pathogenic'
    if any('uncertain' in v for v in vals): return 'VUS'
    if any('benign' in v for v in vals): return 'benign'
    return 'other'

def choose_best_consequence(consequences: List[Dict[str, Any]], canonical_id: Optional[str] = None, gene_symbol_from_query: Optional[str] = None) -> Optional[Dict[str, Any]]:
    """
    Selects the most relevant transcript consequence from a VEP response.
    """
    if not consequences: return None

    if gene_symbol_from_query:
        targeted_consequences = [
            c for c in consequences if c.get('gene_symbol', '').upper() == gene_symbol_from_query.upper()
        ]
        if targeted_consequences:
            consequences = targeted_consequences

    mane_select = [c for c in consequences if c.get('mane_select')]
    if mane_select: return mane_select[0]
    
    if canonical_id:
        canonical_id_base = canonical_id.split('.')[0]
        for c in consequences:
            if c.get('transcript_id', '').startswith(canonical_id_base): return c
            
    coding_consequences = sorted(
        [c for c in consequences if c.get('biotype') == 'protein_coding' and c.get('cds_end')],
        key=lambda c: c['cds_end'] - c['cds_start'] if c.get('cds_start') else -1, reverse=True
    )
    
    return coding_consequences[0] if coding_consequences else consequences[0]

def extract_exons_from_transcript(transcript: Dict[str, Any]):
    exons_raw = sorted(transcript.get('Exon', []), key=lambda e: e['start'])
    if transcript.get('strand') == -1: exons_raw.reverse()
    
    translation_data = transcript.get('Translation', {})
    cds_start, cds_end = translation_data.get('start'), translation_data.get('end')
    seq_region = transcript.get('seq_region_name')
    
    normalized, coding_exon_count = [], 0
    for i, e in enumerate(exons_raw, 1):
        start, end = e['start'], e['end']
        cds_len_of_exon, is_coding = 0, False
        if cds_start and cds_end:
            overlap_start, overlap_end = max(start, cds_start), min(end, cds_end)
            if overlap_end >= overlap_start:
                cds_len_of_exon = overlap_end - overlap_start + 1
                is_coding = True
                coding_exon_count += 1
        
        normalized.append({
            'total_exon_number': i, 'coding_exon_number': coding_exon_count if is_coding else None,
            'exon_id': e.get('id'), 'start': start, 'end': end,
            'seq_region_name': seq_region, 'cds_length': cds_len_of_exon
        })
    return normalized

# --- ASO Strategy Assessment Logic ---

def check_n1c_registry(gene_symbol: str, original_query: str, formatted_hgvs: str) -> Optional[Dict[str, Any]]:
    """
    Searches the pre-loaded N1C registry DataFrame for a matching variant.
    """
    # Check if the DataFrame was loaded successfully
    if n1c_variants_df is None or n1c_variants_df.empty or not gene_symbol:
        return None

    # Extract the core c. notation from the VEP-formatted HGVS string
    core_hgvs_match = re.search(r'(c\..*)', formatted_hgvs, re.IGNORECASE)
    if not core_hgvs_match:
        return None
    core_hgvs = core_hgvs_match.group(1)

    # Filter the DataFrame for the correct gene (case-insensitive)
    gene_matches = n1c_variants_df[n1c_variants_df['Gene'].str.upper() == gene_symbol.upper()]
    if gene_matches.empty:
        return None
        
    # Search within the gene-specific rows for the variant notation
    for index, row in gene_matches.iterrows():
        # Check if our core HGVS notation is present in the registry's 'Coding DNA change (c.)' field
        registry_c_dot = row.get('Coding DNA change (c.)')
        if registry_c_dot and core_hgvs.lower() in registry_c_dot.lower():
            # Match found! Extract data and return the result.
            status = row.get('Status', 'N/A')
            modality = row.get('Therapeutic Modality', 'N/A')
            therapypublication = row.get('Therapy Publication','N/A')
            n1cid = row.get('ID')
            
            return {
                "classification": "Eligible",
                "reason": (
                    f"<p>A direct match for a variant in <strong>{gene_symbol}</strong> was found in the N=1 Collaborative Projects Registry.</p>"
                    f"<ul>"
                    f"<li><strong>Matched Variant:</strong> {registry_c_dot}</li>"
                    f"<li><strong>Status:</strong> {status}</li>"
                    f"<li><strong>Therapeutic Modality:</strong> {modality}</li>"
                    f"</ul>"
                    f"<a href='https://generegistry.n1collaborative.org/entry.html?id={n1cid}' target='_blank' rel='noopener noreferrer'>Click here to view the N=1 Collaborative registry page.</a>"
                ),
                "link": f"https://generegistry.n1collaborative.org/entry.html?id={n1cid}" 
            }
    return None

def _get_c_notation_from_row(row: pd.Series) -> Optional[str]:
    """Best-effort extraction of c. notation from an assessed variants row."""
    candidate_cols = [
        'Variant (c.)', 'Coding DNA change (c.)', 'c_dot', 'HGVS', 'HGVS.c', 'Variant', 'Variant_c'
    ]
    for col in candidate_cols:
        if col in row:
            val = row.get(col)
            if isinstance(val, str) and 'c.' in val:
                return val
    # As a fallback, try to find any field containing 'c.' pattern
    for col in row.index:
        try:
            val = row[col]
            if isinstance(val, str) and 'c.' in val:
                return val
        except Exception:
            continue
    return None

def check_n1c_assessed_variants(gene_symbol: str, formatted_hgvs: str) -> Optional[Dict[str, Any]]:
    """Checks the N1C assessed variants dataset for a curated match and returns a curated assessment."""
    if n1c_assessed_df is None or n1c_assessed_df.empty or not gene_symbol or not formatted_hgvs:
        return None

    core_hgvs_match = re.search(r'(c\..*)', formatted_hgvs, re.IGNORECASE)
    if not core_hgvs_match:
        return None
    core_hgvs = core_hgvs_match.group(1)

    gene_matches = n1c_assessed_df[n1c_assessed_df['Gene'].str.upper() == gene_symbol.upper()] if 'Gene' in n1c_assessed_df.columns else n1c_assessed_df
    if gene_matches.empty:
        return None

    def _normalize_curated_classification(row: pd.Series) -> Tuple[str, Optional[str]]:
        """Attempts to return (classification, reason) based on row fields."""
        # Known target classes for coloring
        allowed = {"Eligible", "Likely Eligible", "Unlikely Eligible", "Not Eligible", "Unable to Assess"}
        # Candidate columns that may encode eligibility
        class_cols = [
            'Eligibility', 'ASO Eligibility', 'Classification', 'classification', 'eligibility', 'Status'
        ]
        val = None
        for col in class_cols:
            if col in row and isinstance(row[col], str) and row[col].strip():
                val = row[col].strip()
                break
        classification = None
        if val:
            v = val.lower()
            if 'not eligible' in v:
                classification = 'Not Eligible'
            elif 'unlikely' in v:
                classification = 'Unlikely Eligible'
            elif 'likely' in v:
                classification = 'Likely Eligible'
            elif 'eligible' in v:
                classification = 'Eligible'
            elif 'unable' in v:
                classification = 'Unable to Assess'
        if classification is None:
            classification = 'Eligible' if 'curated' in (row.get('Tags','') or '').lower() else 'Unable to Assess'
        # Reason/rationale
        reason_cols = ['Notes', 'Rationale', 'Reason', 'Comment']
        reason = None
        for col in reason_cols:
            if col in row and isinstance(row[col], str) and row[col].strip():
                reason = row[col].strip()
                break
        return classification, reason

    for _, row in gene_matches.iterrows():
        c_not = _get_c_notation_from_row(row)
        if c_not and core_hgvs.lower() in str(c_not).lower():
            # Build link using ID when present; prefer variant_entry.html
            link = None
            variant_id_val = row.get('ID')
            if variant_id_val not in (None, "", float('nan')):
                try:
                    link = f"https://generegistry.n1collaborative.org/variant_entry.html?id={variant_id_val}"
                except Exception:
                    link = None
            if not link and 'Link' in row:
                link = row['Link']

            classification, curated_reason = _normalize_curated_classification(row)
            reason = (
                f"<p>A curated assessment exists in the N=1 Collaborative Assessed Variants database for <strong>{gene_symbol}</strong>.</p>"
                f"<p><strong>Matched Variant:</strong> {c_not}</p>"
            )
            if curated_reason:
                reason += f"<p><strong>Notes:</strong> {curated_reason}</p>"
            if link:
                reason += f"<p><a href='{link}' target='_blank' rel='noopener noreferrer'>View N1C curated entry</a></p>"
            return {
                "classification": classification,
                "reason": reason,
                "link": link
            }
    return None

def _extract_exon_numbers_from_text(text: str) -> List[int]:
    """Extracts exon numbers (including ranges) from free text."""
    if not isinstance(text, str) or not text:
        return []
    exons: List[int] = []
    try:
        # Ranges like "exons 45-55" or "exon 2-3"
        for m in re.finditer(r"exons?\s*(\d{1,3})\s*[-–to]+\s*(\d{1,3})", text, flags=re.IGNORECASE):
            a, b = int(m.group(1)), int(m.group(2))
            if a <= b:
                exons.extend(list(range(a, b + 1)))
            else:
                exons.extend(list(range(b, a + 1)))
        # Singles like "exon 51" or "exon-51"
        for m in re.finditer(r"exons?[-\s]*(\d{1,3})", text, flags=re.IGNORECASE):
            try:
                val = int(m.group(1))
                exons.append(val)
            except Exception:
                pass
    except Exception:
        return []
    # Deduplicate while preserving order
    seen = set()
    out: List[int] = []
    for e in exons:
        if e not in seen:
            seen.add(e)
            out.append(e)
    return out

def n1c_exon_skipping_exon_numbers_for_gene(gene_symbol: str) -> Tuple[set, List[str]]:
    """
    Returns (set_of_exon_numbers, list_of_links) for N1C registry rows that indicate exon skipping
    for the given gene. Best-effort extraction across free-text columns.
    """
    exon_set: set = set()
    links: List[str] = []
    if n1c_variants_df is None or n1c_variants_df.empty or not gene_symbol:
        return exon_set, links
    try:
        gene_matches = n1c_variants_df[n1c_variants_df.get('Gene', '').str.upper() == gene_symbol.upper()]
    except Exception:
        # Fallback: if 'Gene' missing due to schema drift
        gene_matches = n1c_variants_df
    if gene_matches is None or len(gene_matches) == 0:
        return exon_set, links

    # Iterate and look for exon skipping hints and extract exon numbers
    for _, row in gene_matches.iterrows():
        row_texts = []
        try:
            for col in row.index:
                val = row[col]
                if isinstance(val, str):
                    row_texts.append(val)
        except Exception:
            pass

        joined = " | ".join(row_texts)
        # Heuristic: must mention exon + skip to avoid false positives
        if re.search(r"exon", joined, re.IGNORECASE) and re.search(r"skip", joined, re.IGNORECASE):
            nums: List[int] = []
            for t in row_texts:
                nums.extend(_extract_exon_numbers_from_text(t))
            for n in nums:
                exon_set.add(n)
            # Build link using ID if available
            try:
                nid = row.get('ID')
                if isinstance(nid, (str, int)) and str(nid).strip() != '':
                    links.append(f"https://generegistry.n1collaborative.org/entry.html?id={nid}")
            except Exception:
                pass

    return exon_set, links

def get_gene_characteristics(gene_symbol: str) -> Dict[str, Any]:
    """
    Retrieves MOI, Haploinsufficiency, and MOA from loaded dataframes,
    including source URLs for better explainability.
    """
    characteristics = {
        "moi": [],
        "haploinsufficiency": {"text": "Unknown", "url": None},
        "moa": [],
        "gene_url": None # --- ADDITION: Initialize gene_url ---
    }

    # --- ClinGen Lookup ---
    if clingen_df is not None and gene_symbol in clingen_df.index:
        # Use .loc[gene_symbol] which can return a Series or DataFrame
        gene_data_row = clingen_df.loc[gene_symbol]
        if isinstance(gene_data_row, pd.DataFrame):
            gene_data_row = gene_data_row.iloc[0] # Take the first row if multiple exist

        gene_url_val = gene_data_row.get('gene_url')
        if pd.notna(gene_url_val):
            characteristics['gene_url'] = gene_url_val

        # 1. Mode of Inheritance (MOI)
        moi_val = gene_data_row.get('mode_of_inheritance')
        if pd.notna(moi_val):
            all_mois = [moi.strip() for moi in str(moi_val).split(',')]
            characteristics["moi"] = sorted(list(set(all_mois)))

        # 2. Haploinsufficiency with Link
        hap_assertion = gene_data_row.get('dosage_haploinsufficiency_assertion')
        hap_url = gene_data_row.get('dosage_report')
        if pd.notna(hap_assertion):
            hap_score_str = str(hap_assertion).strip()
            haplo_text = "Unknown"
            if hap_score_str.startswith('3 -'): haplo_text = "Sufficient evidence"
            elif hap_score_str.startswith('1 -'): haplo_text = "Little evidence"
            elif hap_score_str.startswith('30 -'): haplo_text = "Gene associated with autosomal recessive phenotype"
            else: haplo_text = "No evidence"
            
            characteristics['haploinsufficiency'] = {
                "text": haplo_text,
                "url": hap_url if pd.notna(hap_url) else None
            }

    # --- GOF/LOF Lookup ---
    if goflof_df is not None and gene_symbol in goflof_df.index:
        gene_data = goflof_df.loc[[gene_symbol]]
        if not gene_data[gene_data['LABEL'].str.contains('GOF', na=False)].empty:
            characteristics["moa"].append("GoF")
        if not gene_data[gene_data['LABEL'].str.contains('LOF', na=False)].empty:
            characteristics["moa"].append("LoF")
        if characteristics["moa"]:
            characteristics["moa"] = sorted(list(set(characteristics["moa"])))

    return characteristics

def assess_knockdown(gene_characteristics: Dict[str, Any]) -> Dict[str, Any]:
    """Assesses eligibility for a knockdown strategy."""
    haplo_obj = gene_characteristics.get("haploinsufficiency", {"text": "Unknown"})
    haplo_status_text = haplo_obj.get("text", "Unknown")
    
    moi = gene_characteristics.get("moi", [])
    is_ad_lof = any("Autosomal Dominant" in m for m in moi) and ('LoF' in gene_characteristics.get("moa", []))
    
    checks = {"Gene is not haploinsufficient": haplo_status_text in ["No evidence", "Little evidence"]}
    
    # Determine if haploinsufficiency status is unknown/missing for warning purposes
    haplo_unknown = (haplo_status_text is None) or (str(haplo_status_text).strip().lower() in ["", "unknown", "n/a"])
    
    if is_ad_lof or haplo_status_text == "Sufficient evidence":
        return {
            "classification": "Not Eligible", # Changed to be more definitive
            "reason": "Gene is associated with haploinsufficiency or Autosomal Dominant LoF, therefore knockdown could lead to unintended consequences.",
            "checks": checks
        }
    else:
        # Base reason
        reason = "Gene is not known to be sensitive to haploinsufficiency."
        
        # If haploinsufficiency evidence is unknown/missing, append a warning
        if haplo_unknown:
            reason += " Warning: Haploinsufficiency status is unknown. Therefore knockdown could lead to unintended consequences."
        return {
            "classification": "Likely Eligible",
            "reason": reason,
            "checks": checks
        }
    
def get_overlapping_genes(self, gene_id):
        """Fetches all genes that overlap with a given Ensembl Gene ID."""
        data = self._get(f"/overlap/id/{gene_id}", params={"feature": "gene"})
        return data if isinstance(data, list) else []

def assess_wt_upregulation(client: EnsemblClient, gene_id: str, gene_symbol: str) -> Dict[str, Any]:
    """
    Assesses for WT upregulation by checking for overlapping NATs and by
    searching for a conventional antisense gene name ([GENE]-AS1).
    Trusts the '-AS1' naming convention without a strict biotype check.
    """
    if not gene_id or not gene_symbol:
        return {"classification": "Unable to Assess", "reason": "Missing Gene ID or Symbol."}

    found_antisense_genes = {} # Use a dict to store unique NATs by their ID

    # Collect curated gene-level features from N1C supplementary table
    supp_details: Dict[str, str] = {}
    try:
        if 'n1c_supp_df' in globals() and n1c_supp_df is not None and 'Gene' in n1c_supp_df.columns:
            matches = n1c_supp_df[n1c_supp_df['Gene'].astype(str).str.strip() == str(gene_symbol).strip()]
            def _norm(val: Any) -> str:
                v = str(val).strip().upper()
                if v in ("Y", "YES"): return "Available"
                if v in ("N", "NO"): return "Not available"
                return "Unknown"
            if not matches.empty:
                row = matches.iloc[0]
                supp_details["uORF"] = _norm(row.get('uORF', 'N/A'))
                supp_details["NAT (curated)"] = _norm(row.get('NAT', 'N/A'))
                supp_details["Poison exon (PE)"] = _norm(row.get('PE', 'N/A'))
            else:
                supp_details["uORF"] = "Unknown"
                supp_details["NAT (curated)"] = "Unknown"
                supp_details["Poison exon (PE)"] = "Unknown"
        else:
            supp_details["uORF"] = "Unknown"
            supp_details["NAT (curated)"] = "Unknown"
            supp_details["Poison exon (PE)"] = "Unknown"
    except Exception:
        # Fail-quiet: do not block assessment if supplemental table can't be parsed
        supp_details["uORF"] = "Unknown"
        supp_details["NAT (curated)"] = "Unknown"
        supp_details["Poison exon (PE)"] = "Unknown"

    try:
        # --- Method 1: Search by genomic coordinate overlap ---
        overlapping_genes = client.get_overlapping_genes(gene_id)
        for gene in overlapping_genes:
            if gene.get('biotype') == 'antisense' and gene.get('id') != gene_id:
                found_antisense_genes[gene['id']] = gene

        # --- Method 2: Search by conventional name ([GENE_SYMBOL]-AS1) ---
        antisense_symbol = f"{gene_symbol}-AS1"
        as_gene = client.lookup_symbol(antisense_symbol)
        
        if as_gene:
            found_antisense_genes[as_gene['id']] = as_gene

        # --- Evaluate curated evidence (uORF / NAT / PE) and Ensembl NAT search ---
        has_uorf = supp_details.get("uORF") == "Available"
        has_nat_cur = supp_details.get("NAT (curated)") == "Available"
        has_pe = supp_details.get("Poison exon (PE)") == "Available"
        curated_features = [name for ok, name in [
            (has_uorf, "uORF"),
            (has_nat_cur, "NAT"),
            (has_pe, "Poison exon")
        ] if ok]

        nat_list = []
        nat_names = ""
        nat_ids: List[str] = []
        if found_antisense_genes:
            nat_list = list(found_antisense_genes.values())
            nat_names = ", ".join([nat.get('external_name', nat['id']) for nat in nat_list])
            nat_ids = [nat['id'] for nat in nat_list]

        # Start details with curated statuses
        details = dict(supp_details)
        # Add Ensembl NAT links if any were found
        if nat_list:
            for nat in nat_list:
                nat_name = nat.get('external_name', nat['id'])
                ensembl_link = f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={nat['id']}"
                details[nat_name] = ensembl_link

        if curated_features:
            reason = f"Curated evidence indicates available: {', '.join(curated_features)} for {gene_symbol}."
            if nat_names:
                reason += f" Additionally, Ensembl found antisense transcript(s): {nat_names}."
            return {
                "classification": "Likely Eligible",
                "reason": reason,
                "details": details,
                "antisense_gene_ids": nat_ids,
                "checks": {"Overlapping antisense transcript found": bool(nat_list)}
            }
        else:
            reason = "No curated evidence for uORF, NAT or poison exon availability."
            if nat_names:
                reason += f" However, Ensembl reports antisense transcript(s): {nat_names}."
            return {
                "classification": "Unlikely Eligible",
                "reason": reason,
                "details": details,
                "antisense_gene_ids": nat_ids,
                "checks": {"Overlapping antisense transcript found": bool(nat_list)}
            }

    except Exception as e:
        return {
            "classification": "Unable to Assess",
            "reason": f"An error occurred while searching for antisense transcripts: {e}",
            "checks": {}
        }

def assess_splice_switching(variant_hgvs: str, vep_data: Dict[str, Any], gene_symbol: str) -> Optional[Dict[str, Any]]:
    """
    Assesses a variant for splice-switching potential, adding method and DOI link.
    If not found in the DB, it returns a prompt for user validation.
    """
    if splicevar_df is None or not variant_hgvs or not gene_symbol:
        return None 

    clean_gene = gene_symbol.strip().upper()
    full_canonical_hgvs = variant_hgvs.strip()
    core_hgvs_match = re.search(r'(c\..*)', full_canonical_hgvs, re.IGNORECASE)
    if not core_hgvs_match: return None
    core_canonical_hgvs = core_hgvs_match.group(1).lower()

    gene_rows = splicevar_df[splicevar_df['gene'].str.strip().str.contains(clean_gene, case=False, na=False)]
    if gene_rows.empty:
        # Not found in SpliceVarDB for this gene � check SSCVDB fallback
        if sscvdb_df is not None and not sscvdb_df.empty and vep_data:
            variant_key = _format_sscvdb_variant_id_from_vep(vep_data)
            if variant_key and 'Variant ID' in sscvdb_df.columns:
                if not sscvdb_df[sscvdb_df['Variant ID'].str.strip().str.lower() == variant_key.strip().lower()].empty:
                    details = {
                        "Source Database": "SSCVDB",
                        "Evidence": "Splice-altering reported in SSCVDB",
                        "SSCVDB Gene Page": f"https://sscvdb.io/gene/{gene_symbol}"
                    }
                    return _evaluate_splice_variant_position(variant_hgvs, vep_data, details)
        # --- Return prompt if not in either database ---
        return {
            "classification": "Not in Database",
            "reason": "This variant was not found in SpliceVarDB or SSCVDB. If there is experimental validation (by qPCR or RNA-seq) that this variant is splice-altering, please confirm below.",
            "user_validation_prompt": True
        }

    for index, row in gene_rows.iterrows():
        db_hgvs = str(row.get('hgvs', '')).strip()
        db_hgvs_lower = db_hgvs.lower()
        
        if db_hgvs_lower.endswith(core_canonical_hgvs):
            splice_info = row
            method = str(splice_info.get('method', 'N/A')).strip()
            classification = str(splice_info.get('classification', 'N/A')).strip().lower()
            has_valid_method = method.lower() in ['rna-seq', 'rt-pcr', 'minigene']
            is_splice_altering = classification == 'splice-altering'

            if not (has_valid_method and is_splice_altering):
                continue

            details = {"Confirmation Method": method}
            doi = str(splice_info.get('doi', '')).strip()
            if doi and doi.lower() not in ['na', 'n/a', '']:
                details["Publication"] = f"https://doi.org/{doi}"

            return _evaluate_splice_variant_position(variant_hgvs, vep_data, details)

    # --- If no exact HGVS match was found in SpliceVarDB, try SSCVDB before prompting ---
    if sscvdb_df is not None and not sscvdb_df.empty and vep_data:
        variant_key = _format_sscvdb_variant_id_from_vep(vep_data)
        if variant_key and 'Variant ID' in sscvdb_df.columns:
            if not sscvdb_df[sscvdb_df['Variant ID'].str.strip().str.lower() == variant_key.strip().lower()].empty:
                details = {
                    "Source Database": "SSCVDB",
                    "Evidence": "Splice-altering reported in SSCVDB",
                    "SSCVDB Gene Page": f"https://sscvdb.io/gene/{gene_symbol}"
                }
                return _evaluate_splice_variant_position(variant_hgvs, vep_data, details)

    # --- Return prompt object if still not found ---
    return {
        "classification": "Not in Database",
        "reason": "This variant was not found in SpliceVarDB. If there is experimental validation (by qPCR or RNA-seq) that this variant is splice-altering, please confirm below.",
        "user_validation_prompt": True
    }

def assess_single_exon(client, original_query, transcript, all_exons, target_exon, vep_entry: Dict[str, Any], refseq_id_for_viewer: Optional[str] = None):
    # --- Step 1: Data Gathering and Calculations ---
    gene_id = transcript.get('Parent') 
    transcript_id = transcript.get('id')
    protein_id = transcript.get("Translation", {}).get("id")
    cds_seq = client.get_cds_sequence(transcript_id)
    
    coding_exons = [e for e in all_exons if e['cds_length'] > 0]
    total_coding_exons = len(coding_exons)
    total_cds_len = sum(e['cds_length'] for e in coding_exons)

    if not target_exon.get('coding_exon_number'):
        return {"classification": "Unable to Assess", "reason": f"The variant maps to exon {target_exon['total_exon_number']}, which is non-coding."}
    
    chrom, start, end = target_exon['seq_region_name'], target_exon['start'], target_exon['end']
    variants_in_region = client.overlap_region_variation(chrom, start, end)
    
    clinvar_url = f"https://www.ncbi.nlm.nih.gov/clinvar/?term=GRCh38%3A{chrom}%3A{start}-{end}"

    counts = {'missense': 0, 'inframe_del': 0, 'splice': 0, 'nonsense': 0, 'frameshift': 0, 'benign_splice':0}
    for v in variants_in_region:
        clclass = classify_variant_clinsig(v.get('clinical_significance'))
        conseq = (v.get("consequence_type") or "").lower()
        if clclass == "pathogenic":
            if "missense" in conseq: counts['missense'] += 1
            elif "inframe_deletion" in conseq: counts['inframe_del'] += 1
            elif "splice_donor" in conseq or "splice_acceptor" in conseq: counts['splice'] += 1
            elif "stop_gained" in conseq: counts['nonsense'] += 1
            elif "frameshift" in conseq: counts['frameshift'] += 1
        elif clclass == 'benign' and ("splice_donor" in conseq or "splice_acceptor" in conseq):
            counts['benign_splice'] += 1
    
    exon_cds_len = target_exon['cds_length']
    coding_exon_number = target_exon['coding_exon_number']
    
    # --- Step 2: Condition Checks ---
    cond1_inframe = (exon_cds_len % 3 == 0)
    cond2_no_stop = False
    if cds_seq:
        try:
            cds_map, current_pos = {}, 0
            sorted_coding_exons = sorted(coding_exons, key=lambda x: x['coding_exon_number'])
            for exon in sorted_coding_exons:
                cds_map[exon['coding_exon_number']] = cds_seq[current_pos : current_pos + exon['cds_length']]
                current_pos += exon['cds_length']
            skipped_cds = "".join(cds_map[i] for i in sorted(cds_map.keys()) if i != coding_exon_number)
            if skipped_cds:
                prot = str(Seq(skipped_cds).translate(to_stop=False))
                cond2_no_stop = "*" not in prot[:-1]
        except Exception: cond2_no_stop = False

    cond3_not_terminal = (coding_exon_number is not None and coding_exon_number not in (1, total_coding_exons))
    cond4_small = (exon_cds_len / total_cds_len) < 0.1 if total_cds_len > 0 else False
    
    domains = client.get_domains(protein_id) if protein_id else []
    overlapping_domain_names = []
    if domains:
        cds_pos_start = sum(e['cds_length'] for e in sorted(coding_exons, key=lambda x: x['coding_exon_number']) if e['coding_exon_number'] < coding_exon_number)
        exon_aa_start, exon_aa_end = (cds_pos_start // 3) + 1, ((cds_pos_start + exon_cds_len -1) // 3) + 1
        for d in domains:
            if d.get('start', 0) <= exon_aa_end and d.get('end', 0) >= exon_aa_start:
                overlapping_domain_names.append(d.get('description', d.get('id', 'Unknown Domain')))

    domain_count = len(overlapping_domain_names)
    cond5_no_domain = domain_count == 0
    cond6_missense = counts['missense'] < 3 + counts['frameshift'] + counts['nonsense']
    cond7_splice = counts['splice'] == 0
    cond8_no_inframe_del = counts['inframe_del'] == 0
    cond9_benign_splice = counts['benign_splice'] > 0
    
    # --- Step 3: Classification Logic Chain ---
    classification, reason = "Undetermined", ""
    if cond9_benign_splice:
        classification, reason = "Eligible", "Exon contains benign splice variants, suggesting it may be safely skipped."
    elif not cond3_not_terminal:
        classification, reason = "Not Eligible", "Exon is the first or last coding exon."
    elif not cond1_inframe:
        classification, reason = "Not Eligible", "Exon is out-of-frame, which would disrupt the reading frame."
    elif not cond2_no_stop:
        classification, reason = "Not Eligible", "Skipping this exon is predicted to create a premature stop codon."
    elif not cond7_splice:
        classification, reason = "Not Eligible", f"Exon contains {counts['splice']} pathogenic splice variant(s), indicating exon loss is pathogenic."
    elif not cond4_small and domain_count > 1:
        classification, reason = "Not Eligible", f"Exon overlaps with {domain_count} protein domains and constitutes >10% of the coding region."
    elif not cond8_no_inframe_del:
        classification, reason = "Not Eligible", f"Exon contains {counts['inframe_del']} pathogenic in-frame deletion(s)."
    elif not cond6_missense:
        classification, reason = "Unlikely Eligible", f"Exon is a mutational hotspot with {counts['missense']} pathogenic missense variants."
    elif not cond4_small:
        classification, reason = "Unlikely Eligible", "Exon constitutes >=10% of the protein, risking major functional loss."
    elif not cond5_no_domain:
        classification, reason = "Unlikely Eligible", f"Exon overlaps with {domain_count} protein domain(s)."
    else:
        classification, reason = "Likely Eligible", "Exon meets the primary criteria for a skippable exon."

    # --- Step 4: Visualization Data Generation ---
    visualization_data = None
    try:
        v_chrom, v_start, v_end = vep_entry.get('seq_region_name'), vep_entry.get('start'), vep_entry.get('end')
        if not all([v_chrom, v_start, v_end]): 
            raise ValueError("Missing variant coordinates for visualization.")
        domain_features = []
        if protein_id and domains:
            cds_map = []
            cumulative_cds_len = 0
            is_reverse_strand = transcript.get('strand') == -1
            for exon in sorted(coding_exons, key=lambda x: x['coding_exon_number']):
                cds_len_of_exon = exon['cds_length']
                cds_map.append({
                    'chr': exon['seq_region_name'], 
                    'genomic_start': exon['start'], 
                    'genomic_end': exon['end'], 
                    'transcript_cds_start': cumulative_cds_len + 1, 
                    'transcript_cds_end': cumulative_cds_len + cds_len_of_exon
                })
                cumulative_cds_len += cds_len_of_exon
            
            for domain in domains:
                domain_cds_start, domain_cds_end = (domain['start'] - 1) * 3 + 1, domain['end'] * 3
                for exon_map_entry in cds_map:
                    overlap_start = max(domain_cds_start, exon_map_entry['transcript_cds_start'])
                    overlap_end = min(domain_cds_end, exon_map_entry['transcript_cds_end'])
                    
                    if overlap_start <= overlap_end:
                        offset_start = overlap_start - exon_map_entry['transcript_cds_start']
                        offset_end = overlap_end - exon_map_entry['transcript_cds_start']
                        
                        if not is_reverse_strand:
                            feat_start = exon_map_entry['genomic_start'] + offset_start
                            feat_end = exon_map_entry['genomic_start'] + offset_end
                        else: # On reverse strand, offsets are from the end
                            feat_start = exon_map_entry['genomic_end'] - offset_end
                            feat_end = exon_map_entry['genomic_end'] - offset_start
                            
                        domain_features.append({
                            "chr": exon_map_entry['chr'], 
                            "start": feat_start - 1, 
                            "end": feat_end, 
                            "name": domain.get('description', domain.get('id', 'Domain'))
                        })

        padding = 1000
        visualization_data = {
            "locus": f"{v_chrom}:{max(1, v_start - padding)}-{v_end + padding}",
            "variantTrack": {"name": "Variant", "features": [{"chr": v_chrom, "start": v_start - 1, "end": v_end, "name": vep_entry.get('id', 'Variant')}]},
            "domainTrack": {"name": "Protein Domains", "features": domain_features} if domain_features else None
        }
    except Exception as e:
        import traceback; traceback.print_exc()
        visualization_data = None
        
    # --- FINAL RETURN STATEMENT ---
    return {
        "classification": classification,
        "reason": reason,
        "frac_cds": f"{(exon_cds_len / total_cds_len * 100):.2f}%" if total_cds_len > 0 else "N/A",
        "pathogenic_variant_counts": counts, 
        "domain_count": domain_count,
        "domain_names": overlapping_domain_names,
        "coding_exon_number": coding_exon_number,
        "total_exon_number": target_exon['total_exon_number'], 
        "gene_id": gene_id,
        "transcript_id": transcript_id,
        "clinvar_url": clinvar_url,
        "checks": {
            "Benign splice variant found": cond9_benign_splice, "Is In-Frame": cond1_inframe, 
            "No New Stop Codon": cond2_no_stop, "Not First/Last Exon": cond3_not_terminal,
            "No Pathogenic Splice Variants": cond7_splice, "No Pathogenic In-Frame Deletions": cond8_no_inframe_del,
            "No Domain Overlap": cond5_no_domain, "Low Missense Count": cond6_missense,
            "Is <10% of Protein": cond4_small
        },
        "visualization": visualization_data
    }

def process_single_variant(query: str, client: EnsemblClient, splice_user_input: Optional[str] = None, moa_user_input: Optional[str] = None) -> Dict[str, Any]:
    """
    Contains the complete assessment logic for a single variant query.
    This version is more robust and handles potential unpacking errors.
    """
    try:
        # --- 1. VEP and Consequence Selection ---
        parsed_output = parse_hgvs_query(query)
        if not isinstance(parsed_output, tuple) or len(parsed_output) != 2:
            return {"classification": "Error", "reason": f"Could not parse the input query: '{query}'. Please check the format."}
            
        hgvs_query, gene_symbol_from_query = parsed_output
        if not hgvs_query:
            return {"classification": "Error", "reason": "Invalid input format. Please use a recognized HGVS format (e.g., 'GENE c.123A>G')."}

        vep_data = client.vep_hgvs(hgvs_query)
        if not vep_data or not isinstance(vep_data, list):
            return {"classification": "Unable to Assess", "reason": f"VEP analysis failed for '{hgvs_query}'. The variant may be invalid or not found."}
        
        vep_entry = vep_data[0]
        all_consequences = vep_entry.get('transcript_consequences', [])
        target_consequence = choose_best_consequence(all_consequences, gene_symbol_from_query=gene_symbol_from_query)
        
        if not target_consequence:
            reason = "VEP did not return a consequence."
            if gene_symbol_from_query:
                reason += f" No valid consequence was found for the specified gene '{gene_symbol_from_query}'."
            return {"classification": "Unable to Assess", "reason": reason}

        gene_symbol = target_consequence['gene_symbol']
        definitive_transcript_id = target_consequence['transcript_id']
        gene_id = target_consequence.get('gene_id')

        # --- 2. Get RefSeq ID (for viewer) ---
        refseq_id_for_viewer = None
        mane_consequence = next((c for c in all_consequences if c.get('mane_select')), None)
        if mane_consequence:
            refseq_match = re.search(r'(NM_[0-9]+\.[0-9]+)', mane_consequence['mane_select'])
            if refseq_match:
                refseq_id_for_viewer = refseq_match.group(1)
        if not refseq_id_for_viewer:
            for c in all_consequences:
                if c.get('transcript_id', '').startswith('NM_'):
                    refseq_id_for_viewer = c['transcript_id']
                    break

        # --- 3. Initialize Result & Gene Characteristics ---
        gene_characteristics = get_gene_characteristics(gene_symbol)
        final_result = {
            "summary": {"gene": gene_symbol, "transcript_id": definitive_transcript_id, **gene_characteristics},
            "assessments": {}
        }
        
        # --- 4. N1C Assessed Variants (curated) Check (Exit early if matched) ---
        assessed_match = check_n1c_assessed_variants(gene_symbol, hgvs_query)
        if assessed_match:
            final_result["assessments"]["N1C_Assessed_Variants"] = assessed_match
            return final_result

        # --- 4b. N1C Registry Check (Exit early if matched) ---
        n1c_result = check_n1c_registry(gene_symbol, query, hgvs_query)
        if n1c_result:
            final_result["assessments"]["N1C_Registry_Check"] = n1c_result
            return final_result

        # --- 5. Gene-Level Strategies (Knockdown, WT Upregulation) ---
        # Determine resolved MoA: user selection takes precedence; otherwise a single known MoA
        original_moa_list = list(gene_characteristics.get("moa", []))
        resolved_moa = moa_user_input if moa_user_input in ("GoF", "LoF") else (original_moa_list[0] if len(original_moa_list) == 1 else None)
        # Expose resolved_moa for frontend rendering while preserving original Known Mechanism
        final_result["summary"]["resolved_moa"] = resolved_moa
        if resolved_moa in ("GoF", "LoF"):
            moi = set(gene_characteristics.get("moi", []))
            if resolved_moa == "GoF":
                # Show only Knockdown assessment
                final_result["assessments"]["Allele_Specific_Knockdown"] = assess_knockdown(gene_characteristics)
            elif resolved_moa == "LoF":
                # Show WT Upregulation only when LoF and Autosomal Dominant MOI
                is_lof_ad = any("Autosomal Dominant" in m for m in moi)
                if is_lof_ad:
                    final_result["assessments"]["WT_Upregulation"] = assess_wt_upregulation(client, gene_id, gene_symbol)

        # --- 6. Variant-Specific Strategies (Splice & Exon Skipping) ---
        
        # Define variant type
        consequence_terms = set(target_consequence.get('consequence_terms', []))
        exonic_terms = {'missense_variant', 'stop_gained', 'frameshift_variant', 'synonymous_variant', 'inframe_deletion', 'inframe_insertion','splice_donor_variant', 'splice_acceptor_variant'}
        splice_terms = {'splice_region_variant', }
        is_exonic = any(term in consequence_terms for term in exonic_terms)
        is_splice_region = any(term in consequence_terms for term in splice_terms)
        
        # Run Splice Switching Assessment
        variant_identifier_from_vep = vep_entry.get('input')
        if variant_identifier_from_vep:
            splice_assessment = None
            if splice_user_input == 'yes':
                details = {"Confirmation Method": "User-provided validation (qPCR/RNA-seq)"}
                splice_assessment = _evaluate_splice_variant_position(variant_identifier_from_vep, vep_entry, details)
            elif splice_user_input == 'no':
                splice_assessment = {"classification": "Not Eligible", "reason": "User confirmed no known splice-altering effect."}
            else:
                splice_assessment = assess_splice_switching(variant_identifier_from_vep, vep_entry, gene_symbol)
            
            if splice_assessment:
                final_result["assessments"]["Splice_Switching"] = splice_assessment

        # Run Exon Skipping Assessment *if* variant is exonic/splice
        if is_exonic:
            exon_skip_assessment_added = False
            transcript_data = client.lookup_id_expand(definitive_transcript_id)
            if transcript_data:
                all_exons = extract_exons_from_transcript(transcript_data)
                v_start, v_end = vep_entry['start'], vep_entry['end']
                target_exon = next((ex for ex in all_exons if ex['seq_region_name'] == vep_entry['seq_region_name'] and max(v_start, ex['start']) <= min(v_end, ex['end'])), None)
                
                if target_exon:
                    exon_skip_result = assess_single_exon(client, query, transcript_data, all_exons, target_exon, vep_entry, refseq_id_for_viewer)
                    if "visualization" in exon_skip_result and exon_skip_result["visualization"]:
                        final_result["visualization"] = exon_skip_result.pop("visualization")
                    # N1C registry exon-skipping support: if N1C lists exon skipping for this exon, mark eligible
                    try:
                        n1c_exons, n1c_links = n1c_exon_skipping_exon_numbers_for_gene(gene_symbol)
                    except Exception:
                        n1c_exons, n1c_links = set(), []
                    if target_exon.get('total_exon_number') in n1c_exons:
                        exon_skip_result = dict(exon_skip_result)
                        exon_num = target_exon.get('total_exon_number')
                        note = f"N1C registry lists exon skipping therapy targeting exon {exon_num} in {gene_symbol}. Variants in this exon are considered skippable."
                        prev_reason = exon_skip_result.get('reason', '')
                        exon_skip_result['reason'] = (prev_reason + ' ' + note).strip()
                        exon_skip_result['classification'] = 'Eligible'
                        det = exon_skip_result.get('details', {}) if isinstance(exon_skip_result.get('details'), dict) else {}
                        if n1c_links:
                            det['N1C Exon Skipping Entries'] = ', '.join(n1c_links)
                        det['Supported Exon'] = str(exon_num)
                        exon_skip_result['details'] = det
                    final_result["assessments"]["Exon_Skipping"] = exon_skip_result
                    exon_skip_assessment_added = True
            
            if not exon_skip_assessment_added:
                # Add this block if the logic above fails to add an assessment
                final_result["assessments"]["Exon_Skipping"] = {
                    "classification": "Unable to Assess",
                    "reason": "Variant is exonic or in a splice region, but the target exon could not be determined (e.g., VEP/Ensembl data issue)."
                }

        # 7. Final Fallback
        # This will now only trigger for non-exonic, non-splice, non-gene-strategy variants.
        if not final_result["assessments"]:
            final_result["assessments"]["General_Assessment"] = { 
                "classification": "Unable to Assess", 
                "reason": "Could not determine a primary ASO strategy. The variant is not exonic and no gene-level strategies (Knockdown, WT Upregulation) are applicable." 
            }
        
        return final_result

    except Exception as e:
        import traceback; traceback.print_exc()
        return {"classification": "Error", "reason": f"An unexpected server error occurred: {str(e)}"}

# --- Main Flask Routes ---
app = Flask(__name__)

load_databases()

@app.route('/')
def index(): return render_template('index.html', title="Tool")
@app.route('/about')
def about(): return render_template('about.html', title="About/Methods")
@app.route('/cite')
def cite(): return render_template('cite.html', title="How to Cite")

@app.route('/api_docs')
def api_docs():
    """Serves the API documentation page."""
    return render_template('api_docs.html', title="API Documentation")

@app.route('/api/v1/assess', methods=['GET'])
def api_assess():
    """
    Handles a single variant assessment via a GET request for programmatic access.
    Returns the full assessment data as JSON.
    """
    query = request.args.get('query')
    if not query:
        return jsonify({"error": "The 'query' parameter is required."}), 400

    client = EnsemblClient()
    result = process_single_variant(query, client)
    
    # Provide more specific HTTP status codes based on the outcome
    classification = result.get("classification")
    if classification == "Error":
        return jsonify({"error": result.get("reason", "An internal server error occurred.")}), 500
    if classification == "Unable to Assess":
         return jsonify({"error": result.get("reason", "Could not assess the provided variant.")}), 404
    
    return jsonify(result)

@app.route('/assess', methods=['POST'])
def assess():
    """Handles a single variant assessment request from the frontend."""
    data = request.get_json()
    if not data or 'query' not in data:
        return jsonify({"classification": "Error", "reason": "No query provided."}), 400
    query = data['query']
    splice_input = data.get('splice_user_input', None)
    moa_input = data.get('moa_user_input', None)
    client = EnsemblClient()
    result = process_single_variant(query, client, splice_user_input=splice_input, moa_user_input=moa_input)
    
    return jsonify(result)
@app.route('/batch_assess', methods=['POST'])
def batch_assess():
    if 'file' not in request.files:
        return jsonify({"error": "No file part"}), 400
    file = request.files['file']
    if file.filename == '':
        return jsonify({"error": "No selected file"}), 400

    try:
        if file.filename.endswith('.xlsx'):
            df = pd.read_excel(file, header=None)
        else: # Handles .csv and .txt
            df = pd.read_csv(file, header=None)
    except Exception as e:
        return jsonify({"error": f"Error reading file: {e}"}), 400

    variants = df[0].dropna().astype(str).tolist()
    client = EnsemblClient()
    output_rows = []

    for variant in variants:
        result = process_single_variant(variant, client)
        
        row = {"Variant": variant}
        summary = result.get("summary", {})
        assessments = result.get("assessments", {})

        row["Gene"] = summary.get("gene", "N/A")
        row["MOI"] = ', '.join(summary.get("moi", []))
        row["MOA"] = ', '.join(summary.get("moa", []))
        # Ensembl Transcript ID and link
        transcript_id = summary.get("transcript_id")
        row["Ensembl Transcript"] = transcript_id or "N/A"
        row["Ensembl Transcript Link"] = (
            f"https://www.ensembl.org/Homo_sapiens/Transcript/Summary?t={transcript_id}" if transcript_id else "N/A"
        )
        
        haplo_info = summary.get("haploinsufficiency", {})
        row["Haploinsufficiency"] = haplo_info.get("text", "N/A")
        row["ClinGen Link"] = haplo_info.get("url", "N/A")
        # Duplicate explicitly as curation link for clarity
        row["ClinGen Curation Link"] = haplo_info.get("url", "N/A")
        
        # --- START: NEW COLUMN LOGIC ---

        # 1. Assess if an ASO exists and add the N1C link(s)
        n1c_registry = assessments.get("N1C_Registry_Check", {}) or {}
        n1c_assessed = assessments.get("N1C_Assessed_Variants", {}) or {}
        if n1c_registry or n1c_assessed:
            row["Existing ASO (N1C)"] = "Yes"
            row["N1C Registry Link"] = n1c_registry.get("link", "N/A")
            row["N1C Assessed (Curated) Link"] = n1c_assessed.get("link", "N/A")
        else:
            row["Existing ASO (N1C)"] = "No"
            row["N1C Registry Link"] = "N/A"
            row["N1C Assessed (Curated) Link"] = "N/A"
            
        # 2. Get the Antisense Transcript ID for WT Upregulation
        wt_up = assessments.get("WT_Upregulation", {})
        antisense_ids = wt_up.get("antisense_gene_ids", [])
        row["Antisense Transcript ID"] = ", ".join(antisense_ids) if antisense_ids else "N/A"

        # --- END: NEW COLUMN LOGIC ---

        # Exon Skipping
        skip = assessments.get("Exon_Skipping", {})
        row["Exon Skipping Assessment"] = skip.get("classification", "NA")
        for check, status in skip.get("checks", {}).items():
            row[f"ES Check: {check}"] = status
        # Ensembl exon view and domains (if available)
        if skip:
            gid = skip.get("gene_id")
            tid = skip.get("transcript_id")
            if gid and tid:
                row["Ensembl Exon View Link"] = f"https://www.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g={gid};t={tid}"
            else:
                row["Ensembl Exon View Link"] = "N/A"
            domain_names = skip.get("domain_names") or []
            row["Domains"] = ", ".join(domain_names) if domain_names else "N/A"

        # Splice Correction
        splice = assessments.get("Splice_Switching", {})
        row["Splice Correction Assessment"] = splice.get("classification", "Unable to Assess")
        row["Splicing Validation DOI"] = splice.get("details", {}).get("Publication DOI", "NA")
        # Splicing DB/Source links if available
        splice_details = splice.get("details", {}) if isinstance(splice.get("details", {}), dict) else {}
        splicing_db_link = splice_details.get("SSCVDB Gene Page") or splice_details.get("Publication") or "N/A"
        row["Splicing DB Link"] = splicing_db_link

        # WT Upregulation and Knockdown (assessments remain)
        row["WT-Upregulation"] = wt_up.get("classification", "NA")
        row["Knockdown"] = assessments.get("Allele_Specific_Knockdown", {}).get("classification", "NA")

        # Manual validation needs and dual MoA assessment when MoA unresolved
        manual_needs = []
        summary_moa_list = summary.get("moa", []) or []
        resolved_moa = summary.get("resolved_moa")
        # Splice manual validation prompt
        if splice.get("user_validation_prompt") or (splice.get("classification") in ("Not in Database", "Unable to Assess")):
            manual_needs.append("Splice validation (Variant was not found in SpliceVarDB/SSCVDB and therefore requires user confirmation)")
        # MoA unclear -> assess both and warn
        if not resolved_moa and (len(summary_moa_list) != 1):
            try:
                gof_res = process_single_variant(variant, client, moa_user_input="GoF")
                lof_res = process_single_variant(variant, client, moa_user_input="LoF")
                kd_gof = (gof_res.get("assessments", {}).get("Allele_Specific_Knockdown", {}) or {}).get("classification", "N/A")
                wt_lof = (lof_res.get("assessments", {}).get("WT_Upregulation", {}) or {}).get("classification", "N/A")
                row["Knockdown (GoF)"] = kd_gof
                row["WT-Upregulation (LoF)"] = wt_lof
                row["MoA Dual Assessment Note"] = "Assessed both: use Knockdown if GoF; use WT-Upregulation if LoF."
            except Exception:
                row["Knockdown (GoF)"] = row.get("Knockdown", "N/A")
                row["WT-Upregulation (LoF)"] = row.get("WT-Upregulation", "N/A")
                row["MoA Dual Assessment Note"] = "MoA dual assessment unavailable."
            manual_needs.append("Mechanism selection (GoF vs LoF)")

        row["Manual Validations Needed"] = "; ".join(manual_needs) if manual_needs else "None"

        # Overall Eligibility: highest across all assessment classifications
        def _normalize_class(c: Optional[str]) -> str:
            if not c:
                return "Unable to Assess"
            s = str(c).strip().lower().replace('-', ' ')
            if 'not eligible' in s:
                return 'Not Eligible'
            if 'likely eligible' in s:
                return 'Likely Eligible'
            if 'unlikely eligible' in s:
                return 'Unlikely Eligible'
            if 'eligible' in s:
                return 'Eligible'
            if 'unable to assess' in s or 'not in database' in s:
                return 'Unable to Assess'
            return 'Unable to Assess'

        rank_order = {
            'Not Eligible': 1,
            'Unable to Assess': 2,
            'Unlikely Eligible': 3,
            'Likely Eligible': 4,
            'Eligible': 5,
        }
        best_label = 'Unable to Assess'
        best_score = 0
        for akey, aval in assessments.items():
            if not isinstance(aval, dict):
                continue
            label = _normalize_class(aval.get('classification'))
            score = rank_order.get(label, 2)
            if score > best_score:
                best_score = score
                best_label = label
        row["Overall Eligibility"] = best_label
        
        output_rows.append(row)

    # --- Create and send the Excel file (no changes needed below this line) ---
    if not output_rows:
        return jsonify({"error": "No variants found in file."}), 400
        
    output_df = pd.DataFrame(output_rows)
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        output_df.to_excel(writer, index=False, sheet_name='AVEC_Batch_Results')
    output.seek(0)
    
    return send_file(
        output,
        as_attachment=True,
        download_name='avec_batch_results.xlsx',
        mimetype='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
    )



