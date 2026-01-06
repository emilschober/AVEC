import os
from flask import Flask, request, jsonify, render_template, send_file
import time
import re
import requests
from functools import lru_cache
from decimal import Decimal, getcontext
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
        .pubmed-footer { margin-top: 1em; text-align: center; }
        .pubmed-footer a { text-decoration: none; padding:.4em .75em; border-radius:6px; border:1px solid var(--border); color:#0c63bd; background:#fff }
        .pubmed-footer a:hover { background:#f8fafc }
        /* Domain list styling */
        .domain-list { margin: .35em 0 0 0; padding-left: 1.1em; }
        .domain-list li { margin: 0 0 .2em 0; }
        /* PubMed accordion spacing */
        .pubmed-accordion { margin-top: 1em; }
        .pubmed-accordion .result-header, .pubmed-accordion .result-header h5 { color: #fff; }

        /* Info icon for linking to methods */
        .info-link { margin-left: 8px; color: inherit; text-decoration: none; font-size: 0.8em; opacity: 0.8; }
        .info-link:hover { opacity: 1; }

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
        /* Override checklist for interactive mode */
        .checklist.interactive li { cursor: pointer; user-select: none; transition: background-color 0.1s ease, border-color 0.1s ease; }
        .checklist.interactive li:hover { background-color: #e9ecef; }
        body.dark .checklist.interactive li:hover { background-color: #1e293b; }
        .override-controls { text-align: center; }

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
            // Handle hash for tab navigation on page load
            if (window.location.hash) {
                openTab(window.location.hash.substring(1));
            }
            btn && btn.addEventListener('click', () => {
                const next = document.body.classList.contains('dark') ? 'light' : 'dark';
                localStorage.setItem(key, next);
                apply(next);
            });
        })();
    </script>
    <main class="container">
        {% block content %}{% endblock %}
        <p class="disclaimer">Disclaimer: This tool does not serve as Medical Advice. Although the content of these webpages may partly be provided by individuals in the medical profession, you acknowledge that provision of such content does not create a medical professional-patient relationship and does not constitute an opinion, medical advice, professional diagnosis, service  or treatment of any condition. The access to general information is provided for educational purposes only, through this site and links to other sites. The content is not recommended or endorsed by any doctor or healthcare provider. The information and content provided are not substitutes for medical or professional care, and the information is not to be used in place of a visit, call, consultation or the advice of  your physician or other healthcare  provider. You are liable or responsible for any advice, course of treatment, diagnosis or any other  information, services or product obtained through this site.</p>
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
    .splice-prompt-buttons button.neutral-btn { background-color: #6c757d; }
    .splice-prompt-buttons button.neutral-btn:hover { background-color: #5a6268; }

    .warning { color: #dc3545; font-weight: 600; }
</style>

<h3>AVEC: Automated Variant Eligibility Calculator</h3>
<p>Enter a variant to assess its eligibility for ASO therapy.</p>
<p class="warning">Results depend on underlying data sources and  external services, and may be incomplete or unavailable. Any use of this tool warrants discussion with a biomedical specialist and clinical judgement by a trained physician.</p>
<form id="assessment-form"> <label for="query">Variant:</label> <input id="query" required placeholder="e.g., NM_015427.4:c.1054G>A"> <button type="submit">Assess</button> </form>
<div id="loader">Assessing...</div>
<div id="results"></div>

<hr style="margin: 2em 0;">

<h4 id="batch-toggle" style="cursor: pointer; user-select: none;">
    Batch Processing &#9662;
</h4>

<div id="batch-content" style="display: none;">
    <p>Upload a .csv, .txt, or .xlsx file with one variant per line in the first column.</p>
    <p><a href="/download_batch_template">Download Batch Template (.xlsx)</a></p>
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
// --- Session Management Helpers ---
function saveSession(data, query, spliceInput, moaInput, commentary, esOverrides, kdOverrides, wtOverrides) {
    if (data) sessionStorage.setItem('avec_data', JSON.stringify(data));
    if (query) sessionStorage.setItem('avec_query', query);
    if (spliceInput !== undefined) sessionStorage.setItem('avec_splice', spliceInput || '');
    if (moaInput !== undefined) sessionStorage.setItem('avec_moa', moaInput || '');
    if (commentary !== undefined) sessionStorage.setItem('avec_commentary', commentary || '');
    if (esOverrides !== undefined) sessionStorage.setItem('avec_es_overrides', JSON.stringify(esOverrides));
    if (kdOverrides !== undefined) sessionStorage.setItem('avec_kd_overrides', JSON.stringify(kdOverrides));
    if (wtOverrides !== undefined) sessionStorage.setItem('avec_wt_overrides', JSON.stringify(wtOverrides));
}

function clearSessionUserInputs() {
    sessionStorage.removeItem('avec_splice');
    sessionStorage.removeItem('avec_moa');
    sessionStorage.removeItem('avec_commentary');
    sessionStorage.removeItem('avec_es_overrides');
    sessionStorage.removeItem('avec_kd_overrides');
    sessionStorage.removeItem('avec_wt_overrides');
    sessionStorage.removeItem('avec_data'); // Clear old data on new submit
}

// Restore session on load
window.addEventListener('DOMContentLoaded', () => {
    const savedQuery = sessionStorage.getItem('avec_query');
    const savedData = sessionStorage.getItem('avec_data');
    const savedSplice = sessionStorage.getItem('avec_splice');
    const savedMoa = sessionStorage.getItem('avec_moa');
    const savedCommentary = sessionStorage.getItem('avec_commentary');
    const savedEsOverrides = sessionStorage.getItem('avec_es_overrides');
    const savedKdOverrides = sessionStorage.getItem('avec_kd_overrides');
    const savedWtOverrides = sessionStorage.getItem('avec_wt_overrides');

    if (savedQuery) {
        document.getElementById('query').value = savedQuery;
    }
    
    // Restore global state
    if (savedSplice) window.userSpliceInput = savedSplice;
    if (savedMoa) window.userMoaInput = savedMoa;
    if (savedCommentary) window.userCommentary = savedCommentary;
    if (savedEsOverrides) window.userEsOverrides = JSON.parse(savedEsOverrides);
    if (savedKdOverrides) window.userKdOverrides = JSON.parse(savedKdOverrides);
    if (savedWtOverrides) window.userWtOverrides = JSON.parse(savedWtOverrides);

    if (savedData) {
        try {
            const data = JSON.parse(savedData);
            displayResults(data);
        } catch (e) {
            console.error("Failed to restore session data", e);
        }
    }
});

// Persist user selections across reassessments within the page lifecycle
window.userSpliceInput = null; // 'yes' | 'no' | null
window.userCommentary = '';
window.userMoaInput = null;    // 'GoF' | 'LoF' | null
window.userEsOverrides = {};   // For Exon Skipping
window.userKdOverrides = {};   // For Knockdown
window.userWtOverrides = {};   // For WT Upregulation

document.getElementById('assessment-form').addEventListener('submit', async function(e) {
    e.preventDefault();
    const resultsDiv = document.getElementById('results');
    const loader = document.getElementById('loader');
    const queryVal = document.getElementById('query').value;
    
    resultsDiv.style.display = 'none';
    resultsDiv.innerHTML = '';
    loader.style.display = 'block';
    
    // New query: reset stored user inputs and session storage for inputs
    window.userSpliceInput = null;
    window.userMoaInput = null;
    window.userCommentary = '';
    window.userEsOverrides = {};
    window.userKdOverrides = {};
    window.userWtOverrides = {};
    clearSessionUserInputs();
    
    // Initial assessment payload contains only the query
    const payload = { 
        query: queryVal 
        // splice_user_input is omitted, so backend defaults to DB check
    }; 
    
    try {
        const response = await fetch('/assess', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(payload) });
        const data = await response.json();
        
        // Save to session
        saveSession(data, queryVal, null, null, '', {}, {}, {});
        
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
    const queryVal = document.getElementById('query').value;
    
    resultsDiv.style.display = 'none';
    resultsDiv.innerHTML = '';
    loader.style.display = 'block';

    // Create the payload with the new 'splice_user_input' flag
    const payload = { 
        query: queryVal,
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
        
        // Save session state
        saveSession(data, queryVal, spliceInput, window.userMoaInput, 
            document.getElementById('user-commentary')?.value,
            window.userEsOverrides, window.userKdOverrides, window.userWtOverrides
        );

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
    const queryVal = document.getElementById('query').value;
    
    resultsDiv.style.display = 'none';
    resultsDiv.innerHTML = '';
    loader.style.display = 'block';

    const payload = {
        query: queryVal,
        moa_user_input: moaChoice, // 'GoF' or 'LoF'
        splice_user_input: window.userSpliceInput,
        exon_skipping_overrides: window.userEsOverrides,
        knockdown_overrides: window.userKdOverrides,
        wt_upregulation_overrides: window.userWtOverrides
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
        
        // Save session state
        saveSession(data, queryVal, window.userSpliceInput, moaChoice, 
            document.getElementById('user-commentary')?.value,
            window.userEsOverrides, window.userKdOverrides, window.userWtOverrides
        );

        displayResults(data);
    } catch (error) {
        console.error('Fetch Error:', error);
        displayResults({ classification: 'Error', reason: 'Could not connect to the server.' });
    } finally {
        loader.style.display = 'none';
    }
}

// --- NEW HELPER FUNCTION FOR GUIDELINE OVERRIDES ---
async function reassessWithGuidelineOverrides(strategy, overrides) {
    const resultsDiv = document.getElementById('results');
    const loader = document.getElementById('loader');
    const queryVal = document.getElementById('query').value;

    resultsDiv.style.display = 'none';
    resultsDiv.innerHTML = '';
    loader.style.display = 'block';

    // Update global state
    if (strategy === 'Exon_Skipping') window.userEsOverrides = overrides;
    if (strategy === 'Knockdown') window.userKdOverrides = overrides;
    if (strategy === 'WT_Upregulation') window.userWtOverrides = overrides;

    const payload = {
        query: queryVal,
        moa_user_input: window.userMoaInput,
        splice_user_input: window.userSpliceInput,
        exon_skipping_overrides: window.userEsOverrides,
        knockdown_overrides: window.userKdOverrides,
        wt_upregulation_overrides: window.userWtOverrides
    };

    try {
        const response = await fetch('/assess', { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(payload) });
        const data = await response.json();
        saveSession(data, queryVal, window.userSpliceInput, window.userMoaInput, window.userCommentary, window.userEsOverrides, window.userKdOverrides, window.userWtOverrides);
        displayResults(data);
    } catch (error) {
        console.error('Fetch Error:', error);
    } finally {
        loader.style.display = 'none';
    }
}

// --- HELPER FOR EXONIC PATHOGENICITY ---
async function reassessWithExonicPathogenicityInput(pathogenicityChoice) {
    const resultsDiv = document.getElementById('results');
    const loader = document.getElementById('loader');
    const queryVal = document.getElementById('query').value;

    resultsDiv.style.display = 'none';
    resultsDiv.innerHTML = '';
    loader.style.display = 'block';

    const payload = {
        query: queryVal,
        exonic_pathogenic_user_input: pathogenicityChoice, // 'yes' or 'no'
        // Persist other user choices
        splice_user_input: window.userSpliceInput,
        moa_user_input: window.userMoaInput,
        exon_skipping_overrides: window.userEsOverrides,
        knockdown_overrides: window.userKdOverrides,
        wt_upregulation_overrides: window.userWtOverrides
    };

    try {
        const response = await fetch('/assess', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(payload)
        });
        const data = await response.json();

        // Persist this choice
        window.userExonicPathogenicityInput = pathogenicityChoice;

        // Save session and display
        saveSession(data, queryVal, window.userSpliceInput, window.userMoaInput, window.userCommentary, window.userEsOverrides, window.userKdOverrides, window.userWtOverrides);
        displayResults(data);
    } catch (error) {
        console.error('Fetch Error:', error);
        displayResults({ classification: 'Error', reason: 'Could not connect to the server.' });
    } finally {
        loader.style.display = 'none';
    }
}

// --- EVENT LISTENER (for 'results' div) ---
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
    if (e.target && e.target.id === 'moa-choice-dn') {
        e.preventDefault();
        reassessWithMoaInput('DN');
    }
    // --- NEW: Exonic Pathogenicity Buttons ---
    if (e.target && e.target.id === 'exonic-pathogenicity-yes') {
        e.preventDefault();
        reassessWithExonicPathogenicityInput('yes');
    }
    if (e.target && e.target.id === 'exonic-pathogenicity-no') {
        e.preventDefault();
        reassessWithExonicPathogenicityInput('no');
    }
    // --- Download Report Button ---
    if (e.target && e.target.id === 'download-report-btn') {
        e.preventDefault();
        downloadReport();
    }
    // --- Commit Overrides Button ---
    if (e.target && e.target.classList.contains('commit-overrides-btn')) {
        e.preventDefault();
        const strategy = e.target.dataset.strategy;
        const form = document.getElementById(`override-form-${strategy}`);
        if (!form) return;

        const overrides = {};
        form.querySelectorAll('li[data-check-name]').forEach(li => {
            const checkName = li.dataset.checkName;
            const isPassed = li.classList.contains('pass');
            overrides[checkName] = isPassed;
        });

        reassessWithGuidelineOverrides(strategy, overrides);
    }
    // --- NEW: Interactive Guideline Checklist Item ---
    if (e.target && e.target.closest('.checklist.interactive li')) {
        e.preventDefault();
        const li = e.target.closest('li');
        const strategy = li.dataset.strategy;

        // Toggle pass/fail state
        const isCurrentlyPass = li.classList.contains('pass');
        li.classList.toggle('pass', !isCurrentlyPass);
        li.classList.toggle('fail', isCurrentlyPass);

        // Show the commit button
        const commitBtnContainer = document.getElementById(`override-controls-${strategy}`);
        if (commitBtnContainer) {
            commitBtnContainer.style.display = 'block';
        }
    }
});


// --- displayResults ---
function displayResults(data) {
    const resultsDiv = document.getElementById('results');
    resultsDiv.style.display = 'block';
    resultsDiv.innerHTML = '';

    if (!data || data.classification === "Error" || data.classification === "Unable to Assess") {
        const classificationClass = (data.classification || "Error").toLowerCase().replace(/ /g, '-');
        resultsDiv.innerHTML = `<div class="result-header ${classificationClass}"><h4>${data.classification || "Error"}</h4></div><div class="result-block"><p>${data.reason || "An unknown error occurred."}</p></div>`;
        return;
    }

    let html = '';\n\n    // Precompute MoA-related values for use across sections\n    const assessmentsObj = data.assessments || {};\n    const summaryObj = data.summary || {};\n    let moaList = summaryObj.moa || [];\n    let note = summaryObj.note || '';\n    let resolvedMoa = summaryObj.resolved_moa || window.userMoaInput || null;\n    if (!resolvedMoa) { const m = note.match(/user selection:\\s*(GoF|LoF)/i); if (m) resolvedMoa = m[1]; }\n    const hasN1C = !!(assessmentsObj && (assessmentsObj.N1C_Registry_Check || assessmentsObj.N1C_Assessed_Variants));\n    const hasExactN1CMatch = hasN1C;\n    const hasN1CKnockdown = !!(assessmentsObj && assessmentsObj.N1C_Gene_Knockdown);

    if (data.summary) {
        let geneHTML = data.summary.gene || 'N/A';
        if (data.summary.gene_url) {
            geneHTML = `<a href="${data.summary.gene_url}" target="_blank" rel="noopener noreferrer">${geneHTML}</a>`;
        }
        
        let haploHTML = (data.summary.haploinsufficiency && data.summary.haploinsufficiency.text) || 'N/A';
        if (data.summary.haploinsufficiency && data.summary.haploinsufficiency.url) {
            haploHTML = `<a href="${data.summary.haploinsufficiency.url}" target="_blank" rel="noopener noreferrer">${haploHTML}</a>`;
        }

        let triploHTML = (data.summary.triplosensitivity && data.summary.triplosensitivity.text) || 'N/A';
        if (data.summary.triplosensitivity && data.summary.triplosensitivity.url) {
            triploHTML = `<a href="${data.summary.triplosensitivity.url}" target="_blank" rel="noopener noreferrer">${triploHTML}</a>`;
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
        const proteinEffectHTML = data.summary.protein_effect || 'N/A';
        const rcnvObj = data.summary.rcnv || {};
        const rcnvLink = rcnvObj.url || "#";
        // Create links for the scores
        const pHaploHTML = rcnvObj.pHaplo && rcnvObj.pHaplo !== 'N/A' 
            ? `<a href="${rcnvLink}" target="_blank" rel="noopener noreferrer" title="Collins et al. 2022">${rcnvObj.pHaplo}</a>` 
            : 'N/A';
        const pTriploHTML = rcnvObj.pTriplo && rcnvObj.pTriplo !== 'N/A' 
            ? `<a href="${rcnvLink}" target="_blank" rel="noopener noreferrer" title="Collins et al. 2022">${rcnvObj.pTriplo}</a>` 
            : 'N/A';

        html += `<div class="summary-block"><h4>Query Summary</h4><ul>
                                <li><strong>Gene:</strong> ${geneHTML}</li>
                                <li><strong>Transcript:</strong> ${transcriptHTML}</li>
                                <li><strong>Protein Effect:</strong> ${proteinEffectHTML}</li>
                                <li><strong>Mode of Inheritance:</strong> ${moiHTML}</li>
                                <li><strong>ClinGen Haploinsufficiency:</strong> ${haploHTML}</li>
                                <li><strong>ClinGen Triplosensitivity (ClinGen):</strong> ${triploHTML}</li>
                                <li><strong>pHaplo:</strong> ${pHaploHTML}</li>
                                <li><strong>pTriplo:</strong> ${pTriploHTML}</li> 
                                <li><strong>Orphanet:</strong> ${orphaLinkHTML}</li>
                                <li><strong>ClinGen Pathomechanism:</strong> ${(data.summary.moa && data.summary.moa.join(', ')) || 'N/A'}</li>
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
        const showMoaBlock = (moaList.length === 0 || moaList.length > 1) && !hasExactN1CMatch;
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
        if (moaList.length === 0 && !resolvedMoa && !hasN1C && !hasExactN1CMatch) {
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
        if (!resolvedMoa && !hasExactN1CMatch) {
            html += `<div class="strategy-block">`;
            html += `<div id="toggle-moa" class="result-header unable-to-assess accordion-toggle" style="cursor: pointer; user-select: none;">
                          <h4>Mechanism of Action: Selection Required <span class="chevron">&#9662;</span></h4>
                      </div>`;
            html += `<div id="content-moa" class="result-block" style="display: block;">
                        <p style="font-weight: bold;">There is no single known mechanism of action available. Please select the mechanism relevant to your variant:</p>
                        <div class="splice-prompt-buttons">
                            <button id="moa-choice-gof" type="button" class="neutral-btn">GoF</button>
                            <button id="moa-choice-lof" type="button" class="neutral-btn">LoF</button>
                            <button id="moa-choice-dn" type="button" class="neutral-btn">DN</button>
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
                html += `\n                    <div class=\"strategy-block pubmed-accordion\">\n                        <div id=\"toggle-moa-pubmed\" class=\"result-header unable-to-assess accordion-toggle\" style=\"cursor: pointer; user-select: none;\">\n                            <h5>PubMed: Mechanism of Action <span class=\"chevron\">&#9662;</span></h5>\n                        </div>\n                        <div id=\"content-moa-pubmed\" class=\"result-block\" style=\"display: none;\">\n                            <p><strong>Related PubMed Search:</strong> ${moaQuery}</p>\n                            <div id=\"pubmed-moa-results\" class=\"pubmed-context\" data-query=\"${moaQuery}\">Loading PubMed results...</div>\n                            <div class=\"pubmed-footer\"><a href=\"${pmUrl}\" target=\"_blank\" rel=\"noopener noreferrer\">Open full results on PubMed</a></div>\n                            <p class=\"note\" style=\"margin-top:0.5em;\">Powered by NCBI E-utilities.</p>\n                        </div>\n                    </div>`;
            })();
            html += `</div></div>`;
        }
        for (const [strategy, result] of Object.entries(data.assessments)) {
            if (hasN1CKnockdown && strategy === 'Knockdown') continue;
            const strategyName = strategy.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase());

            // --- START: MODIFIED LOGIC FOR SPLICE_SWITCHING & EXONIC PATHOGENICITY (COMBINED) ---
            if (strategy === 'Splice_Switching' && result.user_validation_prompt_exonic_pathogenicity) {
                // Render the special "Awaiting Pathogenicity Input" block
                html += `<div class="strategy-block">`;
                html += `<div id="toggle-${strategy}" class="result-header unable-to-assess accordion-toggle" style="cursor: pointer; user-select: none;">
                             <h4>${strategyName}: Awaiting User Input <span class="chevron">&#9662;</span></h4>
                         </div>`;
                // Set to 'block' so the user sees the question immediately
                html += `<div id="content-${strategy}" class="result-block" style="display: block;">`; 
                html += `<p style="font-weight: bold;">${result.reason}</p>`;
                // The new buttons
                html += `<div class="splice-prompt-buttons">
                             <button id="exonic-pathogenicity-yes" type="button" class="no-btn">Pathogenic Effects (independent of splicing)</button>
                             <button id="exonic-pathogenicity-no" type="button" class="yes-btn">No Pathogenic Effects (independent of splicing)</button>
                         </div>`;
                html += `</div></div>`;
            
            } else if (strategy === 'Splice_Switching' && result.user_validation_prompt) {
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
                    const terms = '(Splicing OR Splice-altering OR "intron retention" OR "splice-switch" OR "cryptic splice site activation")';
                    const spliceQuery = prefix ? `${prefix} AND ${terms}` : terms;
                    html += `<div class="pubmed-footer" style="margin-top:1em;"><a href="https://pubmed.ncbi.nlm.nih.gov/?term=${encodeURIComponent(spliceQuery)}" target="_blank" rel="noopener noreferrer">Search PubMed for Splicing Evidence</a></div>`;
                })();
                html += `</div></div>`;
            } else {
                // This is the "else" block: it contains all the *original* rendering logic
                let classificationClass = (result.classification || "Unable to Assess").toLowerCase().replace(/ /g, '-');
                // --- START: CUSTOM CLASSIFICATION MAPPING ---
                const strategyToTabMap = {
                    'Exon_Skipping': 'skipping',
                    'Splice_Switching': 'splice',
                    'Knockdown': 'knockdown',
                    'WT_Upregulation': 'wt'
                };
                const tabId = strategyToTabMap[strategy];
                const infoLink = tabId ? `<a href="/about#${tabId}" class="info-link" title="View methods for ${strategyName}">&#9432;</a>` : '';
                // Map new WT Upregulation labels to existing CSS classes for color consistency
                if (classificationClass === 'potential-possibilities-identified') {
                    classificationClass = 'likely-eligible';
                } else if (classificationClass === 'no-potential-possibilities-identified') {
                    classificationClass = 'unlikely-eligible';
                } else if (classificationClass === 'manual-assessment-required') {
                    // NEW: Map 'Manual Assessment Required' to the 'Unable to Assess' style
                    classificationClass = 'unable-to-assess';
                }
                // --- END: CUSTOM CLASSIFICATION MAPPING ---

                
                html += `<div class="strategy-block">`;
                // DUPLICATE HEADER REMOVED HERE
                html += `<div id="toggle-${strategy}" class="result-header ${classificationClass} accordion-toggle" style="cursor: pointer; user-select: none;">
                             <h4>${strategyName} ${infoLink} <span class="chip-status status-${classificationClass}">${result.classification || "N/A"}</span> <span class="chevron">&#9662;</span></h4>
                         </div>`;
                // Set to 'none' by default
                html += `<div id="content-${strategy}" class="result-block" style="display: none;">`;
                
                if (strategy === 'Exon_Skipping' && result.total_exon_number && result.gene_id && result.transcript_id) {
                    const ensemblLink = `https://www.ensembl.org/Homo_sapiens/Transcript/Exons?db=core;g=${result.gene_id};t=${result.transcript_id}`;
                    html += `<p><strong>Target:</strong> <a href="${ensemblLink}" target="_blank" rel="noopener noreferrer">Exon ${result.total_exon_number}</a></p>`;
                }
                if (strategy === 'Exon_Skipping' && result.exon_length) {
                    html += `<p><strong>Exon Length:</strong> ${result.exon_length} nt</p>`;
                }
                if(result.reason) html += `<p><strong>Reason:</strong> ${result.reason}</p>`;
                
                if(result.note) {
                    html += `<p class="warning" style="margin-top: 0.5em;"><strong>Warning:</strong> ${result.note}</p>`;
                }
                
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

                if (strategy === 'Splice_Switching') {
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
                        const spliceToggleId = `toggle-${strategy}-pubmed`;
                        const spliceContentId = `content-${strategy}-pubmed`;
                        const spliceResultsId = `pubmed-${strategy}-results`;
                        html += `\n                        <div class=\"strategy-block pubmed-accordion\">\n                            <div id=\"${spliceToggleId}\" class=\"result-header unable-to-assess accordion-toggle\" style=\"cursor: pointer; user-select: none;\">\n                                <h5>PubMed: Splicing Evidence <span class=\"chevron\">&#9662;</span></h5>\n                            </div>\n                            <div id=\"${spliceContentId}\" class=\"result-block\" style=\"display: none;\">\n                                <p><strong>Related PubMed Search:</strong> ${spliceQuery}</p>\n                                <div id=\"${spliceResultsId}\" class=\"pubmed-context\" data-query=\"${spliceQuery}\">Loading PubMed results...</div>\n                                <div class=\"pubmed-footer\"><a href=\"${pmUrl}\" target=\"_blank\" rel=\"noopener noreferrer\">Open full results on PubMed</a></div>\n                                <p class=\"note\" style=\"margin-top:0.5em;\">Powered by NCBI E-utilities.</p>\n                            </div>\n                        </div>`;
                    })();
                }

                if (result.checks && strategy !== 'Knockdown' && strategy !== 'WT_Upregulation' && (result.classification || '').toLowerCase() !== 'eligible') {
                    const totalChecks = Object.keys(result.checks).length; 
                    const passedChecks = Object.values(result.checks).filter(c => c.passed).length;
                    const hasEffectiveOverrides = Object.values(result.checks).some(c => c.overridden && c.passed !== c.original_passed);
                    let overrideNote = hasEffectiveOverrides ? ' <span class="note">(Manual Override Applied)</span>' : '';
                    
                    html += `<h5>Guideline Checks <span class="chip ${passedChecks===totalChecks ? 'chip-ok':'chip-warn'}">${passedChecks}/${totalChecks} passed</span>${overrideNote}</h5>`;
                    html += `<div id="override-form-${strategy}" class="override-form"><ul class="checklist interactive">`;
                    for (const [check, status] of Object.entries(result.checks)) {
                        const isPassed = status.passed;
                        const isChanged = status.overridden && (status.passed !== status.original_passed);
                        const overrideIndicator = isChanged ? ' <span class="note">(override)</span>' : '';
                        const textContent = `${check}${overrideIndicator}`;
                        html += `<li class="${isPassed ? 'pass' : 'fail'}" data-check-name="${check}" data-strategy="${strategy}">
                                     ${isChanged ? `<strong>${textContent}</strong>` : textContent}
                                 </li>`;
                    }
                    html += '</ul></div>';
                }
                html += `<div id="override-controls-${strategy}" class="override-controls" style="display: none; margin-top: 1em;"><button type="button" class="commit-overrides-btn" data-strategy="${strategy}">Commit Manual Overrides</button></div>`;
                
                if (strategy === 'Knockdown' && result.checks) {
                    const totalChecks = Object.keys(result.checks).length;
                    const passedChecks = Object.values(result.checks).filter(c => c.passed).length;
                    const hasEffectiveOverrides = Object.values(result.checks).some(c => c.overridden && c.passed !== c.original_passed);
                    let overrideNote = hasEffectiveOverrides ? ' <span class="note">(Manual Override Applied)</span>' : '';

                    html += `<h5>Guideline Checks <span class="chip ${passedChecks === totalChecks ? 'chip-ok' : 'chip-warn'}">${passedChecks}/${totalChecks} passed</span>${overrideNote}</h5>`;
                    html += `<div id="override-form-${strategy}" class="override-form"><ul class="checklist interactive">`;
                    for (const [check, status] of Object.entries(result.checks)) {
                        const isPassed = status.passed;
                        const isChanged = status.overridden && (status.passed !== status.original_passed);
                        const overrideIndicator = isChanged ? ' <span class="note">(override)</span>' : '';
                        const textContent = `${check}${overrideIndicator}`;
                        html += `<li class="${isPassed ? 'pass' : 'fail'}" data-check-name="${check}" data-strategy="${strategy}">${isChanged ? `<strong>${textContent}</strong>` : textContent}</li>`;
                    }
                    html += '</ul></div>';
                    html += `<div id="override-controls-${strategy}" class="override-controls" style="display: none; margin-top: 1em;"><button type="button" class="commit-overrides-btn" data-strategy="${strategy}">Commit Manual Overrides</button></div>`;
                }
                if (result.pathogenic_variant_counts) {
                    let evidenceHeader = '<h5>Evidence from Databases</h5>';
                    if (result.clinvar_url) {
                        evidenceHeader = `<h5><a href="${result.clinvar_url}" target="_blank" rel="noopener noreferrer">Evidence from Databases (ClinVar)</a></h5>`;
                    }
                    html += evidenceHeader;

                    let domainHTML = result.domain_count || 'N/A';
                    if (Array.isArray(result.domain_details) && result.domain_details.length > 0) {
                        const domainItems = result.domain_details.map(d => {
                            const baseLabel = d.label || [d.name, d.source].filter(Boolean).join(' ') || 'Domain';
                            const label = d.source && !baseLabel.includes(d.source) ? `${baseLabel} (${d.source})` : baseLabel;
                            if (d.url) return `<li><a href="${d.url}" target="_blank" rel="noopener noreferrer">${label}</a></li>`;
                            return `<li>${label}</li>`;
                        }).join('');
                        domainHTML = `<ul class="domain-list">${domainItems}</ul>`;
                    } else if (result.domain_names && result.domain_names.length > 0) {
                        domainHTML = result.domain_names.join(', ');
                    }

                    const spliceVariantLinks = result.splice_variant_links || {};
                    const renderVariantLinks = (variants) => {
                        if (!Array.isArray(variants) || variants.length === 0) return '';
                        const items = variants.map(v => {
                            const label = v.id || 'Variant';
                            if (v.url) return `<a href="${v.url}" target="_blank" rel="noopener noreferrer">${label}</a>`;
                            return label;
                        }).join(', ');
                        return `<div class="note" style="margin: 0.15em 0 0 0.75em;">${items}</div>`;
                    };
                    const pathogenicSpliceLinks = renderVariantLinks(spliceVariantLinks.pathogenic);
                    const benignSpliceLinks = renderVariantLinks(spliceVariantLinks.benign);
                    const benignSpliceCount = result.pathogenic_variant_counts.benign_splice || 0;

                    html += `<ul>
                                    <li><strong>Fraction of Protein:</strong> ${result.frac_cds || 'N/A'}</li>
                                    <li><strong>Overlapping Protein Domains:</strong> ${domainHTML}</li>
                                    <li><strong>Pathogenic Variants in Exon:</strong><ul>
                                        <li>Missense: ${result.pathogenic_variant_counts.missense}</li>
                                        <li>Nonsense: ${result.pathogenic_variant_counts.nonsense}</li>
                                        <li>Frameshift: ${result.pathogenic_variant_counts.frameshift}</li>
                                        <li>In-frame Deletions: ${result.pathogenic_variant_counts.inframe_del}</li>
                                        <li>Splice Site: ${result.pathogenic_variant_counts.splice}${pathogenicSpliceLinks}</li>
                                        <li>Benign Splice Site: ${benignSpliceCount}${benignSpliceLinks}</li>
                                    </ul></li>
                                </ul>`;
                }

                if (strategy === 'Knockdown') {
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
                        const terms = '((antisense OR ASO) AND (knockdown OR \"gene silencing\" OR \"allele specific\" OR \"allele-specific\"))';
                        const kdQuery = prefix ? `${prefix} AND ${terms}` : terms;
                        const pmUrl = 'https://pubmed.ncbi.nlm.nih.gov/?term=' + encodeURIComponent(kdQuery);
                        const kdToggleId = `toggle-${strategy}-pubmed`;
                        const kdContentId = `content-${strategy}-pubmed`;
                        const kdResultsId = `pubmed-${strategy}-results`;
                        html += `\n                        <div class=\"strategy-block pubmed-accordion\">\n                            <div id=\"${kdToggleId}\" class=\"result-header unable-to-assess accordion-toggle\" style=\"cursor: pointer; user-select: none;\">\n                                <h5>PubMed: Knockdown Evidence <span class=\"chevron\">&#9662;</span></h5>\n                            </div>\n                            <div id=\"${kdContentId}\" class=\"result-block\" style=\"display: none;\">\n                                <p><strong>Related PubMed Search:</strong> ${kdQuery}</p>\n                                <div id=\"${kdResultsId}\" class=\"pubmed-context\" data-query=\"${kdQuery}\">Loading PubMed results...</div>\n                                <div class=\"pubmed-footer\"><a href=\"${pmUrl}\" target=\"_blank\" rel=\"noopener noreferrer\">Open full results on PubMed</a></div>\n                                <p class=\"note\" style=\"margin-top:0.5em;\">Powered by NCBI E-utilities.</p>\n                            </div>\n                        </div>`;
                    })();
                }
                html += `</div></div>`;
            }
            // --- END: MODIFIED LOGIC ---
        }
    }
    
    // Add Commentary & Download section at the end
    html += `<div class="strategy-block">
                <div id="toggle-report" class="result-header unable-to-assess accordion-toggle" style="cursor: pointer; user-select: none;">
                    <h4>Add Commentary & Download Report <span class="chevron">&#9662;</span></h4>
                </div>
                <div id="content-report" class="result-block" style="display: none;">
                    <h5>Your Commentary</h5>
                    <textarea id="user-commentary" placeholder="Add your notes, interpretation, or next steps here. This will be included in the downloaded report." style="width: 100%; min-height: 120px; padding: .5em; border-radius: 4px; border: 1px solid #ccc; font-family: inherit; font-size: .95em;"></textarea>
                    <button id="download-report-btn" type="button" style="margin-top: 1em;">Download HTML Report</button>
                </div>
            </div>`;
    resultsDiv.innerHTML = html;\n\n    // Expand/Collapse all wiring\n    (function(){\n        const expandBtn = document.getElementById('expand-all');\n        const collapseBtn = document.getElementById('collapse-all');\n        if (expandBtn) expandBtn.addEventListener('click', () => {\n            document.querySelectorAll('.result-block').forEach(el => el.style.display='block');\n            document.querySelectorAll('.accordion-toggle').forEach(t => t.classList.add('open'));\n        });\n        if (collapseBtn) collapseBtn.addEventListener('click', () => {\n            document.querySelectorAll('.result-block').forEach(el => el.style.display='none');\n            document.querySelectorAll('.accordion-toggle').forEach(t => t.classList.remove('open'));\n        });\n    })();

    // Restore commentary from session/global state
    const commentaryTextarea = document.getElementById('user-commentary');
    if (commentaryTextarea && window.userCommentary) {
        commentaryTextarea.value = window.userCommentary;
    }
    // Save commentary to global state on input
    if(commentaryTextarea) commentaryTextarea.addEventListener('input', (e) => { window.userCommentary = e.target.value; });

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
                       `<div class="pubmed-meta">${journal}${journal && pubdate ? '  ' : ''}${pubdate}${authors ? '  ' + authors : ''}</div>` +
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
        const targets = document.querySelectorAll('.pubmed-context, #pubmed-moa-results');
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

function downloadReport() {
    const query = document.getElementById('query').value;
    const commentary = document.getElementById('user-commentary').value;
    const timestamp = new Date().toLocaleString();

    let reportHTML = `
        <!doctype html>
        <html lang="en">
        <head>
            <meta charset="utf-8">
            <title>AVEC Report: ${query}</title>
            <style>
                body { font-family: system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif; line-height: 1.6; color: #333; max-width: 800px; margin: 2em auto; padding: 0 1em; }
                h1, h2, h3 { color: #005a9c; border-bottom: 1px solid #eee; padding-bottom: 0.3em; }
                a { color: #005a9c; }
                .summary-block, .assessment-block { border: 1px solid #ddd; border-radius: 8px; padding: 1em; margin-bottom: 1.5em; background-color: #f9f9f9; }
                .summary-block ul { list-style: none; padding: 0; }
                .summary-block li { margin-bottom: 0.5em; }
                .assessment-block h3 { margin-top: 0; font-size: 1.2em; }
                .disclaimer { font-size: 0.8em; color: #777; text-align: center; margin-top: 2em; }
                .commentary { background-color: #fffbeb; border-color: #fde68a; }
                pre { background-color: #f0f0f0; border: 1px solid #ddd; border-radius: 4px; padding: 1em; white-space: pre-wrap; word-wrap: break-word; }
            </style>
        </head>
        <body>
            <h1>AVEC Assessment Report</h1>
            <p><strong>Query:</strong> ${query}</p>
            <p><strong>Report Generated:</strong> ${timestamp}</p>
    `;

    if (commentary) {
        reportHTML += `<div class="summary-block commentary"><h2>User Commentary</h2><pre>${commentary}</pre></div>`;
    }

    const summaryContent = document.querySelector('.summary-block');
    if (summaryContent) {
        reportHTML += `<h2>Query Summary</h2><div class="summary-block">${summaryContent.innerHTML}</div>`;
    }

    reportHTML += `<h2>Therapeutic Assessments</h2>`;
    document.querySelectorAll('.strategy-block').forEach(block => {
        // Exclude the commentary/download block itself from the report
        if (block.querySelector('#toggle-report')) return;
        const header = block.querySelector('.result-header h4');
        const content = block.querySelector('.result-block');
        if (header && content) {
            reportHTML += `<div class="assessment-block"><h3>${header.innerText.replace('▼', '').replace('▲', '')}</h3>${content.innerHTML}</div>`;
        }
    });

    reportHTML += `<p class="disclaimer">This tool is for informational purposes only and is not for clinical use. Results require manual verification.</p></body></html>`;

    const blob = new Blob([reportHTML], { type: 'text/html' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = `AVEC_Report_${query.replace(/[^a-zA-Z0-9]/g, '_')}.html`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(a.href);
}
</script>
{% endblock %}
"""

    about_html = """
{% extends "base.html" %}
{% block content %}
<style>
    /* Styling für das Tab-System */
    .tab-nav {
        display: flex;
        flex-wrap: wrap;
        gap: 0;
        margin-bottom: 20px;
        border-bottom: 1px solid var(--border);
    }
    .tab-btn {
        background: transparent;
        border: none;
        border-bottom: 3px solid transparent;
        color: var(--muted);
        padding: 10px 20px;
        cursor: pointer;
        font-weight: 600;
        font-size: 1em;
        transition: all 0.2s ease;
        border-radius: 4px 4px 0 0;
        box-shadow: none; /* Override default button shadow */
        width: auto; /* Override default full width */
        margin: 0;
    }
    .tab-btn:hover {
        color: var(--text);
        background-color: rgba(0,0,0,0.03);
    }
    .tab-btn.active {
        color: var(--primary);
        border-bottom-color: var(--primary);
        background-color: transparent;
    }
    .tab-pane {
        display: none;
        animation: fadeIn 0.3s ease-in-out;
        line-height: 1.6;
    }
    .tab-pane.active {
        display: block;
    }
    .tab-pane h4 {
        margin-top: 0;
        margin-bottom: 1em;
        color: var(--text);
        border-bottom: 1px solid var(--border);
        padding-bottom: 0.5em;
    }
    @keyframes fadeIn {
        from { opacity: 0; transform: translateY(5px); }
        to { opacity: 1; transform: translateY(0); }
    }
    /* Dark mode adjustments */
    body.dark .tab-btn:hover { background-color: rgba(255,255,255,0.05); color: var(--text); }
</style>

<h3>About & Methods</h3>

<div class="tab-nav">
    <button class="tab-btn active" onclick="openTab('general')">General</button>
    <button class="tab-btn" onclick="openTab('splice')">Splice Correction</button>
    <button class="tab-btn" onclick="openTab('skipping')">Exon Skipping</button>
    <button class="tab-btn" onclick="openTab('knockdown')">Knockdown</button>
    <button class="tab-btn" onclick="openTab('wt')">WT Upregulation</button>
</div>

<div id="general" class="tab-pane active">
    <h4>Purpose & Methodology</h4>
    <p>AVEC (Automated Variant Eligibility Calculator) is a research tool designed to assess the eligibility of a variant for therapeutic antisense oligonucleotides (ASOs). It automates the analysis of key criteria to predict whether ASO therapy (splice-correcting, WT-upregulation, knockdown or exon skipping) is likely to restore a functional protein product, and assumes variants are nonsense, missense small indels or frameshift variants. For other variants please refer to the <a href="https://shorturl.at/YqphL" target="_blank" rel="noopener noreferrer"><b>guidelines</b></a>.</p>
    <p>If the underlying datasources are false or incomplete it is possible that the results are also false. Therefore we strongly advise manual validation of the results.</p>
    <p>This is a prototype.</p>
</div>

<div id="splice" class="tab-pane">
    <h4>Evaluation of Splicing Effects</h4>
    <p>
        In this step, assessors will evaluate the effect on splicing of each variant. As explained in the Background, correction of aberrant splicing is an elegant way to restore the reading frame and the physiological splicing pattern (Fig. 1). Thus, the first aim should always be to assess whether a splice-switching ASO strategy is applicable. Whenever this is not a possibility, the flowchart in Step 4 will aid in choosing the most suitable section for variant assessment. For the splicing evaluation and correction of aberrant splicing, different considerations have to be made given inheritance patterns and pathomechanism of the genetic variant (described separately at the end of Step 3). Table 3 provides an overview of the classification of variants for eligibility for ASO splice correction.
    </p>

    <div class="toc-block" style="background: #f8f9fa; padding: 1em; border-left: 4px solid #005a9c; margin-bottom: 2em;">
        <strong>In this section, the following is covered:</strong>
        <ul>
            <li>Determining whether a variant affects splicing and what is considered sufficient evidence for mis-splicing</li>
            <li>Considerations for eligibility of splice correction ASOs for intronic and exonic (including UTR) variants</li>
            <li>Types of exonic variants that can cause aberrant splicing and alternate ASO strategies to consider in place of splice correction</li>
            <li>Special considerations for canonical splice site variants</li>
            <li>Important considerations for pathomechanism and inheritance pattern</li>
            <li>How to search the literature for splice correction ASOs</li>
        </ul>
    </div>

    <p>
        Different types of variants can influence the splicing process and thus the decision on which ASO strategy to apply. While all variants should be assessed for their splice-altering potential (Anna & Monika, 2018), some exceptions usually do not have to be evaluated in Step 3:
    </p>
    <ul>
        <li>For nonsense and frameshift variants that are associated with a LoF mechanism at the protein level, see Section A (or Section C in cases of haploinsufficiency). Rare cases where nonsense and frameshift variants that affect splicing lead to GoF or DN effects can be assessed with the considerations outlined in Step 3.</li>
        <li>For whole exon duplications and deletions, please consider the variant effect and jump to Section B (knockdown) or Section C (upregulation from the wildtype), check Fig. 7 for directions.</li>
    </ul>

    <p>
        The effects on splicing of each variant need to be confirmed with a functional assay, whereby only RNAseq, qPCR, or cDNA sequencing/analysis obtained from patient-derived cells can be considered sufficiently reliable. Please note that it is considered sufficient if functional data of the above-mentioned kind is available on the same variant from a different individual, i.e., a case report in the literature or information provided in ClinVar. Data gathered through artificial assays, like exon-trapping and mini- and midi-genes, cannot be considered sufficient, as these assays do not recapitulate endogenous splicing, and results can be misleading (Lin et al., 2021). Explicitly, the prediction of splicing effects using in silico tools is not sufficient evidence and cannot be used for these assessments (Oh et al., 2024).
    </p>

    <p>
        In case there is functional evidence against a splice-altering effect, i.e., confirmation that the variant does not affect splicing, please refer to Fig. 7 to guide you towards the next section applicable to your assessment. Additionally, if there is no functional evidence of splicing effects, please refer to Fig. 7. In cases where there is no evidence of splicing effects (whether for or against), this does not mean the variant is ineligible for a splice correcting ASO. A lack of functional evidence towards splicing means this variant cannot be assessed for slice correction eligibility at this time. Theoretically, all variants are suspicious of splicing until demonstrated otherwise, and if evidence on splicing effects becomes available, the variant should be reassessed for eligibility towards splice correction ASOs.
    </p>

    <p>
        If splicing effects are confirmed, the exact effects need to be evaluated. This can be a gain of an acceptor or donor splice site or the loss/weakening of an acceptor or donor splice site. Especially, if canonical splicing is destroyed - due to canonical splice site variants, branchpoint variants, or variants destroying other splice-regulatory elements - the variant is most likely not eligible for a splice correction ASO treatment. Canonical splicing is considered destroyed if no wildtype transcript/splicing at the canonical splice sites can be identified in the functional analysis. If some wildtype transcript is still produced or protein function is detected, canonical splice sites are considered as weakened. Variants that fully abolish canonical splicing are usually not amenable to a splice correction ASO and are classified as “not eligible” for this approach. For other ways to target these variants, please see our section on canonical splice site variants.
    </p>

    <p>
        We recommend studying available data carefully and to also assess data provided in supplementary material, as one can often find relevant gels and blots in the supplement that may only be hinted at in the main manuscript. Further, look out for evidence illustrated by gels or qPCR results etc. instead of solely relying on how the authors describe or discuss the data. From our experience, some manuscripts claim there is no wildtype splicing left whereas faint wildtype bands can be identified in a gel, which might be sufficient to consider this a suitable case. Please consider that a variant's effect on splicing may not have been assessed. The absence of evidence might indicate that splicing effects were overlooked.
    </p>

    <p>
        We generally distinguish intronic and exonic variants for the splicing assessments. Please also see Fig. 5 for a general overview of the classification of variants towards splice correction.
    </p>

    <h5>Intronic variants</h5>
    <p>
        As a rule of thumb, we consider intronic variants that are >-100 bp (upstream of acceptor splice site) or >+50 bp (downstream of donor splice site) away from the nearest canonical splice site as likely eligible, i.e., there should be no negative impact on the canonical splice sites and branchpoint (Fig. 5). For variants closer to the canonical splice sites, enough functional evidence should be available on the effect on canonical splicing. Variants within 5 bp of the exon-intron/intron-exon boundary, often referred to as splice-site variants, are considered “not eligible” for splice correction. There are some instances where splice-site variants lead to leaky exon skipping (Dawes et al., 2023; Feng et al., 2022). If functional data shows that the variant is causing leaky exon skipping, meaning that the exon is only skipped in some RNA transcripts, and the canonical splice site is weakened but not destroyed, please refer to the section titled “Considerations for Exon Inclusion” at this step (Step 3) for further information. Note that the variant should also not disturb the branchpoint. The location of the branchpoint can vary and is usually found 18-40 bp upstream from the acceptor splice site, but can also be found anywhere 10-100 bp upstream (i.e., -10 to -100 bp region) (Xie, Wang, & Lin, 2023).
    </p>
    <p>
        Intronic variants with functional evidence of causing aberrant splicing can now be analyzed using the above-mentioned criteria and classified using Table 3.
    </p>

    <img src="{{ url_for('static', filename='fig5.png') }}" class="responsive-img" alt="Eligible targets for splice correcting ASOs">
    <p class="img-caption">
        <strong>Figure 5: Eligible targets for splice correcting ASOs</strong><br>
        The design of splice correcting ASOs takes into consideration key splicing motifs. The most amenable variants for splice correction are deep intronic variants, highlighted yellow in this figure. As the canonical splice sites are approached, one must consider the effects of the ASO blocking important splice site motifs (i.e., the branchpoint). These regions are highlighted in black (intronic) or gray (exonic) in this figure. Anything within 5 bp is considered not eligible for targeting with an ASO (highlighted red in this figure, both intronic and exonic), though exceptions apply.
    </p>

    <h5>Exonic variants</h5>
    <p>
        For exonic variants, further considerations apply. As for the intronic variants, the exonic variant should be at least 5 bp outside of the canonical splice sites (upstream and downstream - hard cutoff). In some instances, a variant within 5bp of the canonical splice site may lead to leaky exon skipping, where the canonical splice signal is weakened but not destroyed. If functional evidence shows that this is the case, refer to the section titled “Considerations for Exon Inclusion” at this step (Step 3). In cases where the variant is destroying the canonical splice site, leading to full exon skipping, these should be treated the same as a complete exonic deletion. Skipping an adjacent exon to restore the reading frame may be possible. For these cases, refer to Section A (exon skipping). We further define a second cut-off at 15 bp (region 6 to 15 bp upstream or downstream of the canonical splice sites) as a soft cutoff where a splice correction can be considered but is challenging. This cutoff is based on the idea that ASOs can bind to this region without destroying the canonical splice set, yet there is still a possibility of weakening canonical splicing.
    </p>
    <p>
        For exonic variants, the specific type of variant is important for further analysis.
    </p>
    <p><strong>Nonsense and frameshift variants that cause aberrant splicing can be further distinguished:</strong></p>
    <ul>
        <li>If the splice aberration itself leads to a LoF effect on the protein, the variant cannot be analyzed using this part of the guidelines. Correction of splicing would still lead to a LoF on the protein level. Thus, these variants fall under Section A for exon skipping analysis (or Section C in the case of haploinsufficiency).</li>
        <li>If the splice aberration leads to a toxic GoF effect on protein level, the variant can be considered for splice correction. Correction of the splicing effect will then lead to an early truncation and most likely a LoF effect on protein level. For certain cases, this is a useful ASO approach. It is also possible to consider these variants for downregulation (Section B).</li>
        <li>In both situations, our considerations for haploinsufficiency should be taken into account.</li>
    </ul>
    
    <p><strong>Synonymous variants</strong> that cause aberrant splicing do not influence the amino acid sequence and can be assessed using Table 3.</p>

    <p><strong>For missense variants and small in-frame indels (&lt;50 bp) (Mahmoud et al., 2019)</strong>, the considerations are more complex as both the effect on splicing and also the effect of altering the protein coding sequence need to be taken into account. It needs to be established that the variant effect at the protein level solely arises from the aberrant splicing and not that a missense variant or small indel itself causes a pathogenic effect. It is possible that a missense variant causing aberrant splicing leading to a LoF at the protein level, could independently cause a GoF or DN effect as a missense variant.</p>
    <ul>
        <li>If the effect of a missense variant or in-frame indel on the protein sequence is not known (independent of the splice-altering effect), the variant cannot be assessed until further evidence is available.</li>
        <li>If the change in amino acid sequence is known to be pathogenic, for example, if a different pathogenic nucleotide change causes the same amino acid change without the splicing effect, different sub-criteria apply:
            <ul>
                <li>Aberrant splicing and missense variants cause GoF/DN effect → can be considered for downregulation (Section B) or for exon skipping (Section A).</li>
                <li>Aberrant splicing and missense variants cause LoF effect, with the amino acid change causing a complete loss of function on the protein level → can be considered for exon skipping (see Section A).</li>
                <li>Aberrant splicing causes GoF effect and missense variant would lead to a LoF → can be considered for exon skipping (Section A) or knockdown (Section B) and also for splice correction if the LoF phenotype is milder or loss of one allele is tolerated.</li>
                <li>Aberrant splicing causes LoF effect and missense variant causes GoF/DN effect → can be considered for exon skipping (Section A) or for splice correction in case of a GoF effect that leads to a less severe phenotype compared to the LoF phenotype.</li>
                <li>Aberrant splicing leads to GoF or LoF effect and the variant causes partial loss of function/reduced protein function → can be considered for splice correction and would be classified as “unlikely eligible” for splice correction. Here, the underlying thought is that having a partially functional protein is better than having no protein function or toxic protein function due to aberrant splicing. The variant can also be considered for exon skipping (Section A) or for downregulation (Section B) in case of GoF effects from aberrant splicing.</li>
            </ul>
        </li>
        <li>If the amino acid change is known to be benign on protein level, due to population data or functional studies, the variant follows the same considerations as a synonymous variant and can be analyzed using Table 3.</li>
    </ul>

    <h5>UTR Variants</h5>
    <p>
        UTR Variants should be assessed in a manner similar to exonic coding variants. If there is functional evidence that a variant in the 5’ or 3’ UTR leads to aberrant splicing, they could be considered for splice correcting ASOs. For example, a 3’ UTR variant in KLHL40 resulted in the creation of a cryptic donor site and was assumed to cause nonsense mediated decay (Dofash et al., 2022). This variant could be treated as an exonic gain-of-splice variant and could be targeted with a splice-correcting ASO.
    </p>

    <h5>Important considerations for different inheritance patterns and pathomechanisms</h5>
    <p>
        For AR disorders, Step 3 applies without restrictions. If the variants are homozygous or compound heterozygous, the Step 3 evaluation has to be done once or twice for both variants. If the disorder is AD associated with a LoF, the Step 3 guidelines also apply. Since the purpose of these ASOs is to restore wildtype splicing, rather than alter splicing, whether the ASO binds the wildtype or pathogenic allele does not matter (unlike with canonical exon skipping ASOs and gapmer ASOs). Similarly, variants that lead to a GoF or DN effect through altering splicing can be assessed with these guidelines (please consider the effect of the variant on protein level once splicing is corrected). For GoF or DN variants, however, one might consider a knockdown approach for a more efficacious effect (see Section B).
    </p>

    <h5>ASO check</h5>
    <p>
        In addition to the recommended assessment strategies, assessors should review the literature for splice correcting ASO strategies. This review can be performed either as the final step to validate the assessment strategy or earlier in the assessment process (see Step 0). Specifically, for splice correcting ASOs, it is crucial to evaluate whether a splice correcting approach has been implemented and validated for the specific variant at the RNA level. Please also consider that the rescue at the RNA level needs to be sufficient to produce enough protein to rescue the phenotype. Ideally, a publication should provide this evidence.
    </p>
    <p>
        Splicing mechanisms differ by variant, making it crucial to ensure that any existing ASOs found in the literature are applicable to the specific variant of interest. Additionally, if the literature indicates that an ASO is ineffective in correcting splicing for the variant of interest, this does not necessarily render the variant ineligible for ASO development. We recommend that a variant should only be considered “not eligible” for splice-correcting ASOs if there is evidence from two independent investigations at the protein or functional level, or from one investigation providing a convincing explanation as to why an ASO cannot be developed.
    </p>
    <p>
        In cases of conflicting evidence, consider the quality of the research, the types of experiments conducted, the nature of the results shared, and the publication date. The evaluation of literature on existing ASOs should be carried out at the assessor's discretion, with a critical and discerning approach.
    </p>
    <p>
        To help with identifying available ASOs for splice correction, we recommend a search term in Pubmed like this (searching for an ASO targeting a specific variant), text in bold would need to be adjusted for the gene and variant being examined:
        <br><code>ABCA4 AND ((ASO) OR (AON) OR (antisense oligonucleotide)) AND ((Gln876Ter) OR (c.2626C>T) OR (E876X) OR (Q876X) OR (Gln876*) OR (E876*) OR (Q876*))</code>
    </p>

    <h5>Table 3: Classification of variants for their eligibility towards splice correction</h5>
    <style>
        .classification-table { width: 100%; border-collapse: collapse; margin-top: 1em; font-size: 0.95em; }
        .classification-table th, .classification-table td { border: 1px solid #ddd; padding: 8px 12px; vertical-align: top; }
        .classification-table th { background-color: #f1f1f1; text-align: left; width: 25%; }
        .bg-green { background-color: #e8f5e9; }
        .bg-blue { background-color: #e1f5fe; }
        .bg-orange { background-color: #fff3e0; }
        .bg-red { background-color: #ffebee; }
    </style>
    <table class="classification-table">
        <thead>
            <tr>
                <th>Classification</th>
                <th>Criteria</th>
            </tr>
        </thead>
        <tbody>
            <tr class="bg-green">
                <td><strong>Eligible</strong></td>
                <td>ASO has already been developed and shown to work with available functional evidence at the protein level (pre-clinical data is sufficient)</td>
            </tr>
            <tr class="bg-blue">
                <td><strong>Likely eligible</strong></td>
                <td>
                    Functional studies (RNAseq, qPCR) validate alternate splicing<br>
                    <strong>AND (if intronic)</strong><br>
                    {<br>
                    Intronic variant -101 bp or +51 bp outside of the canonical splice sites<br>
                    OR<br>
                    No weakened branch point/canonical splice site as determined by functional studies<br>
                    }<br>
                    <strong>AND (if exonic, including UTR)</strong><br>
                    Donor/acceptor gained is not within 15 bp of a canonical splice site<br>
                    <strong>AND (if exonic coding)</strong><br>
                    Functional evidence shows there is no pathogenic effect of the amino acid change(s) on the protein.
                </td>
            </tr>
            <tr class="bg-orange">
                <td><strong>Unlikely eligible</strong></td>
                <td>
                    Functional studies (RNAseq, qPCR) validate alternate splicing<br>
                    <strong>BUT</strong><br>
                    {<br>
                    Canonical splice site and branchpoint is weakened (but still functional, i.e., either canonical transcript or protein function still detectable)<br>
                    OR/AND (if intronic)<br>
                    Intronic variant within -6 and -100 bp or +6 and +50 bp of the canonical splice sites.<br>
                    OR/AND (if exonic)<br>
                    Donor/acceptor gained within 6-15bp of the canonical splice site.<br>
                    OR/AND (if exonic)<br>
                    There is evidence of residual protein function (i.e., if splicing is corrected, the protein coding change still produces a partially functional protein).<br>
                    }
                </td>
            </tr>
            <tr class="bg-red">
                <td><strong>Not eligible</strong></td>
                <td>
                    Canonical splice site and branchpoints are destroyed - functional evidence shows exon skipping in 100% of transcripts on the mutated allele<br>
                    OR<br>
                    Different nucleotide change leading to the same predicted amino acid residue change – but no alternate splicing – is pathogenic<br>
                    OR<br>
                    Variant within 5 bp of the canonical splice site<br>
                    OR<br>
                    Evidence that ASO cannot be developed, shown by two independent investigations on the protein/functional level or one investigation with a convincing explanation of why ASO cannot be developed.
                </td>
            </tr>
        </tbody>
    </table>

    <hr style="margin: 2em 0; border: 0; border-top: 1px solid #ddd;">

    <h4>Considerations for Exon Inclusion</h4>
    <h5>Introduction to Exon Inclusion</h5>
    <p>
        As described earlier in Step 3, canonical splice site variants and variants disrupting splice regulatory elements (e.g., Branchpoints, splice enhancers/silencer) are excluded from further assessment for splice correction, and are classified as not eligible. However, in some instances, weakened splice sites leading to exon skipping can be strengthened through blocking a splice silencer signal, meaning silencing the silencer will lead to increased exon inclusion. Exon inclusion is the mechanism behind the most widely used ASO to date, nusinersen, that enhances the inclusion of exon 7 into the SMN2 transcript (NM_017411.4) used to treat patients with spinal muscular atrophy (Hua et al., 2008). This is currently also the only drug in clinical application using this strategy, highlighting the difficulty in achieving sufficient exon inclusion. Further examples of successful pre-clinical exon inclusion have been published for SLC26A4, PAH, and GAA (Feng et al., 2022; Martínez-Pizarro et al., 2024; van der Wal et al., 2017).
    </p>
    <p>
        In this section of Step 3, assessors will be given considerations on how to assess towards exon inclusion. This can also include alternative exons that are rarely spliced in or exons that have a lower inclusion frequency due to pathogenic variants weakening the splice sites. The latter case is especially interesting as inclusion of these exons would lead to restoration of the canonical transcripts. At the same time, this is where most often the challenge lies, exons that are skipped due to pathogenic variants are incredibly difficult to include back into the transcripts using ASOs.
        It is not recommended to perform exon inclusion assessment in a routine setting and instead should be performed by assessors skilled in this type of assessment. However, for completeness, we discuss how these assessments are generally performed. Unlike other ASO strategies discussed in these guidelines (splice correction, exon skipping, and knockdown) we do not assess eligibility for these variants here. We do however consider a variant eligible for exon inclusion in cases where an exon inclusion strategy has been published with clear functional data that goes beyond mini- and midi-gene studies, similar to considerations for wildtype upregulation (Section C).
    </p>

    <h5>Availability of Functional Data</h5>
    <p>
        To assess a variant for potential exon inclusion, functional data on the missplicing caused by the variant is necessary. There needs to be a clear understanding of the effect of the variant and the different transcripts produced. Preferably, this data is obtained from patient-derived cells. Please note that a mini-/midi-gene assay alone is not sufficient. Furthermore, the level of wildtype splicing needs to be established. Exon inclusion should only be considered for variants that weaken canonical splicing, either by weakening the canonical splice site or branchpoint directly or by disrupting a splice enhancer. Generally, there should still be substantial wildtype splicing present in the mutant transcripts. From the published literature where exon inclusion was shown preclinically, usually at least 20% (and up to 70%) of wildtype transcript was still present in the patients (Feng et al., 2022; Hua et al., 2008; Martínez-Pizarro et al., 2024; van der Wal et al., 2017). We thus recommend to ensure that sufficient levels of wildtype transcript are still produced.
    </p>
    <p>
        Fully destroyed canonical splice sites or instances where canonical splicing is fully abolished cannot be rescued with exon inclusion strategies. Our explanations described earlier in this section (Step 3) are still true for these considerations here. These variants are considered “not eligible” for splice correction.
    </p>

    <h5>Exon Inclusion ASO Target Sites</h5>
    <p>
        When the exon skipping is caused by disruption of an exonic splice enhancer, exon inclusion can only be promoted by blocking an exonic or intronic splice silencer working in tandem. Similar considerations apply when exon skipping is caused by a canonical splice site variant. Such counteracting motifs are not always present.
    </p>
    <p>
        Most often hnRNP regulatory motifs are being targeted, ideally in the intron. These motifs can be identified through searching for the sequences in the intronic region close to the respective splice site, usually an ASO walk is then performed to identify an ASO sequence that best blocks the silencer. Here, it should be ensured that the motif is not too close to the canonical splice site so that blocking the silencer with an ASO is not negatively impacting splicing on the canonical splice site. It should also be made sure that the splice silencer is indeed affecting the splice site in question and not regulating a different canonical or cryptic splice site. These analyses and studies are beyond the scope of these guidelines.
    </p>
    <p>
        To summarize, for exon inclusion strategies, splicing data must be available from a patient-derived sample (endogenous splicing) that clearly indicates the levels of wildtype and missplicing. Should a significant amount of wildtype transcript still be present, we recommend reaching out to an expert group to inquire for assessment for exon inclusion strategies, including where to target.
    </p>

    <hr style="margin: 2em 0; border: 0; border-top: 1px solid #ddd;">

    <h4>Canonical Splice Site Variants and Splice Region Variants</h4>
    <p>
        As discussed above, canonical splice site variants and variants in close proximity to the canonical splice sites (with 5 bp) are usually not good candidates for an ASO strategy. The reason being that most often the canonical splicing is disturbed making it hard to correct these effects. Given the high demand for ASO therapies involving these pathogenic variants and the fact that exceptions do apply, we here discuss scenarios in which canonical splice site variants can be considered for an ASO approach. We recommend discussing such an approach with an expert before starting the development.
    </p>
    <p>
        Pathogenic variants at or around the canonical splice sites can lead to a multitude of effects including full exon skipping, leaky exon skipping, intron retention and cryptic splicing. It is also possible that a pathogenic canonical splice site variants leads to double- or even multi-exon skipping. An overview of potential effects is given in Fig. 6 and the matching ASO strategy is provided in Table 4. Please note that Fig. 6 is not exhaustive and combinations of effects can occur.
    </p>
    <p>
        Before an ASO strategy can be chosen, functional data of the consequences of the variant need to be obtained that also include the level of wildtype transcript that is still present. Without such data, no proper strategy can be determined. In such cases, the variant remains “unable to assess”. Importantly, it can also be that a canonical splice site variant does not lead to splice disruption or minimal disruption not being disease-causing. Assessors should ensure that the variant is indeed the cause of disease before targeting such a variant.
    </p>

    <img src="{{ url_for('static', filename='fig6.png') }}" class="responsive-img" alt="Possible aberrant splicing effects of canonical splice site variants">
    <p class="img-caption">
        <strong>Figure 6: Possible aberrant splicing effects of canonical splice site variants</strong><br>
        A) Full exon skipping can lead to different effects depending on whether the exon is in-frame our out-of-frame. B) In the case of leaky exon skipping, canonical mRNA transcript is still produced. C) In case of cryptic exon splicing due to a canonical splice site variant, the cryptic splice site can be in- frame our out-of-frame with the canonical splice site leading to different effects. D) Partial intron retention can occur when intronic cryptic splice sites are used instead of the canonical splice site that is disrupted. In some cases (E), the full intron is retained. In cases C-E, the effects can also be leaky with some canonical splicing still being present.
    </p>

    <h5>Table 4: Overview of consequences of canonical splice site variants and potential therapeutic strategies</h5>
    <table class="classification-table">
        <thead>
            <tr>
                <th>Effect on splicing</th>
                <th>In-frame exon</th>
                <th>Out-of-frame exon</th>
            </tr>
        </thead>
        <tbody>
            <tr>
                <td><strong>Full exon skipping</strong></td>
                <td>Not eligible for splice correction or exon skipping approaches, consider upregulation from the wildtype in case of LoF with haploinsufficiency or knockdown if the effect is DN or Go</td>
                <td>Can be treated as a single-exon deletion and adjacents exons can potentially be skipped. Further consider upregulation from the wildtype in case of LoF with haploinsufficiency or knockdown if the effect is DN or GoF</td>
            </tr>
            <tr>
                <td><strong>Partial (leaky) exon skipping</strong></td>
                <td colspan="2">If sufficient WT splicing is left one can consider exon inclusion IF the effect can be counteracted by blocking a strong splice silencer. Further consider upregulation from the wildtype in case of LoF with haploinsufficiency or knockdown if the effect is DN or GoF.</td>
            </tr>
            <tr>
                <td><strong>Cryptic splicing exonic</strong></td>
                <td>Not eligible for splice correction or exon skipping as the loss of an in-frame part of the exon leads to the phenotype. Consider upregulation from the wildtype in case of LoF with haploinsufficiency or knockdown if the effect is DN or GoF. If cryptic splicing is leaky, see leaky exon skipping.</td>
                <td>Can be considered for full exon skipping (Section A) in case where the exon itself is in-frame and fulfills the criteria for exon skipping. Consider upregulation from the wildtype in case of LoF with haploinsufficiency or knockdown if the effect is DN or GoF. If cryptic splicing is leaky, see leaky exon skipping.</td>
            </tr>
            <tr>
                <td><strong>Cryptic splicing intronic/partial intron retention</strong></td>
                <td>Usually not eligible for splice correction if the canonical splice site is fully destroyed. If splicing is leaky, consider blocking the cryptic splice site with an ASO or try to enforce full exon skipping. Consider upregulation from the wildtype in case of LoF with haploinsufficiency or knockdown if the effect is DN or GoF.</td>
                <td>Usually not eligible for splice correction if the canonical splice site is fully destroyed. If splicing is leaky, consider blocking the cryptic splice site with an ASO. Consider upregulation from the wildtype in case of LoF with haploinsufficiency or knockdown if the effect is DN or GoF.</td>
            </tr>
            <tr>
                <td><strong>Full intron retention</strong></td>
                <td colspan="2">Usually not eligible for splice correction and exon skipping. Consider upregulation from the wildtype in case of LoF with haploinsufficiency or knockdown if the effect is DN or GoF.</td>
            </tr>
        </tbody>
    </table>
</div>

<div id="skipping" class="tab-pane">
    <h4>Considerations for Canonical Exon Skipping</h4>
    <p>
        In this section, assessors will be guided on how to assess a variant for eligibility towards canonical exon skipping ASOs. As described in the Background, exon skipping ASOs can be used to skip exons containing pathogenic variants (Fig. 1). Variants generally to be considered for exon skipping are nonsense and frameshift variants causing a LoF in AR, AD, X-linked and Y-linked disorders (additional considerations apply in dominantly inherited diseases, see “Important considerations for different inheritance patterns and pathomechanisms” at the end of this section). GoF and DN negative variants can also be amenable for an exon skipping approach under certain circumstances as well as whole exon deletions (see “Important considerations for different inheritance patterns and pathomechanisms” in this section). While splicing assessments are based on the variant, evaluating the potential for exon skipping are based on the exon, thus an exon skipping approach is applicable to all variants within that exon. For an exon to be skippable, different criteria need to be met (Fig. 8 ).
    </p>

    <div class="toc-block" style="background: #f8f9fa; padding: 1em; border-left: 4px solid #005a9c; margin-bottom: 2em;">
        <strong>This section covers the following topics:</strong>
        <ul>
            <li>Assessing exon position and frame</li>
            <li>Assessing exon size</li>
            <li>Strategies for searching and considering naturally occurring exon skipping and in-frame deletions</li>
            <li>Assessing the role of functional domains</li>
            <li>Important considerations for different inheritance patterns and pathomechanisms</li>
            <li>Strategies for searching the literature for exon skipping ASOs</li>
        </ul>
    </div>

    <img src="{{ url_for('static', filename='fig8.png') }}" class="responsive-img" alt="Overview of exon skipping assessment">
    <p class="img-caption">
        <strong>Figure 8: Overview of exon skipping assessment using a hypothetical transcript.</strong><br>
        The first and last coding exons (exons 2 and 10) cannot be skipped, in addition to out-of-frame exons (exons 3,4,7,9, shown by the shape of the exons). In-frame exons encoding for more than 10% of the coding region (exon 5) and exons coding for functional domains (exon 8) are unlikely eligible for exon skipping. Likely eligible for exon skipping are small, in-frame exons that do not code for a domain (e.g., exon 6). White: exon not considered for canonical exon skipping, e.g., exons in the UTRs, grey: exons not eligible for skipping, blue: exons unlikely eligible for exon skipping, green: exon likely eligible for exon skipping, orange: functional domain.
    </p>

    <h5>Assessing Exon Position and Frame</h5>
    <p>
        The first and last coding exons are usually not eligible for exon skipping as this would lead to the loss of the start or stop codon. Note, due to untranslated regions (UTRs) the first and last exons of a transcript are not necessarily coding exons and, thus, assessors should be aware of which exons contain the canonical start and stop codons (Aspden, Wallace, & Whiffin, 2023) (Fig. 8). Additionally, genes containing exactly one coding exon are not eligible and are disqualified from further analysis. Next, assessors should evaluate the length and position/frame of the exon. Out-of-frame exons (containing a number of base pairs not divisible by 3) should be classified as “not eligible” in the case of LoF variants where restoration of the reading frame is the aim. Assessors can use tools such as the UCSC Genome Browser, ExonViz (van den Berg, Lauffer, & Laros, 2024), or Ensembl to determine whether the exon is in-frame or out-of-frame and whether the exon is coding or non-coding (Fig. 9).
    </p>

    <img src="{{ url_for('static', filename='fig9.png') }}" class="responsive-img" alt="Determining exon frames">
    <p class="img-caption">
        <strong>Figure 9: Determining exon frames.</strong><br>
        A) Exon frames in ExonViz can be identified by their shape. Rectangular shapes are in-frame with a phase at 0-0 (exon starts and ends with the codon). Exons with an arrow on one end and a notch on the other end are also considered in-frame, whereby it does not matter if the arrow is upstream or downstream. If there are two arrows or two notches, the exons are out-of-frame. An arrow indicates a 1 nucleotide overhang while a notch illustrates a 2 nucleotide overhang. Small, in-frame exons are circled in red (exon 6, 7, 17 and 25). For example, exon 6 starts with an arrow, meaning the exon starts with 1 nucleotide from a codon starting in the previous exon (i.e., this is the third nucleotide in the codon, with the first 2 nucleotides being found in exon 5). Exon 6 ends with a notch, indicating that these are the first two nucleotides of a codon (2-2 phase). In total, both ends together cover 3 nucleotides. B) Identification with Ensembl can be done by checking the phases of an exon. Exons with phases 0-0, 1-1, 2-2 are in-frame. Additionally, dividing the exon length by 3 can help determine the frame. C) UCSC Genome Browser shows phases to determine the frames using mouse over. The browser directly indicates the phase for you. Please ensure you are hovering over the correct transcript.
    </p>

    <p>
        For in-frame exons, it is possible that the first and last codon is partially encoded by the adjacent exons (exons with phases 1-1 and 2-2). In this case, the formation of a new codon (including stop codons) is possible at the new exon-exon boundary. Assessors should evaluate whether joining the adjacent exons would form an in-frame stop codon, coded by either a TAA, TAG, or TGA (Fig. 10). If this is the case, the exon is not eligible for exon skipping (exceptions apply, please see below). If a novel codon is created that results in a different amino acid at the the same residue, we recommend classifying this as unlikely eligible. Assessors can utilize the UCSC Genome Browser or Ensembl to check (instructions on how to perform this check are provided in example video 11, the NM_000170.3(GLDC):c.538C>T [MIM 238300] example).
    </p>

    <img src="{{ url_for('static', filename='fig10.png') }}" class="responsive-img" alt="Formation of a new codon">
    <p class="img-caption">
        <strong>Figure 10: Formation of a new codon as a result of exon skipping.</strong><br>
        When an in-frame exon with phase 1-1 or 2-2 is skipped, the adjacent exons will form a new codon on the new exon-exon boundary. A) shows the formation of a stop codon when an ASO is used to skip in-frame exon 2. The “T” nucleotide from exon 1, and “AA” nucleotides from exon 3 join together to form a stop codon. This leads to a premature termination resulting in an absent gene product or a truncated, non-functional product. B) shows the formation of a new codon with a functional product. The “GG” nucleotides at the end of exon 1 join the “A” nucleotide at the start of exon 3, upon skipping of exon 2. This codes for glycine.
    </p>

    <p>
        The skipping of an exon that disrupts the reading frame (e.g., skipping of an out-of-frame exon, formation of a new termination codon when an exon is skipped, or skipping the first and last coding exon) tends to be not eligible for analysis. However, there are exceptions to these rules:
    </p>
    <ul>
        <li>It is theoretically possible to skip out-of-frame exons to generate a premature termination codon within the last or penultimate exon given that there is some left-over function of the protein (please note that this is theoretically possible - one would still have to assess the last exon for the importance of domains and whether this shortened transcript is predicted to undergo nonsense mediated decay)</li>
        <li>It is also possible to skip out-of-frame exons for GoF or DN variants to downregulate transcript levels, please see “Considerations for different inheritance patterns and pathomechanisms” in this section, as well as Section B</li>
        <li>It is theoretically possible to skip an in-frame exon which results in the formation of a stop codon if the resulting protein is stable and has some left-over function (unlikely eligible)</li>
        <li>It is possible to skip adjacent out-of-frame exons in case of the variant being a whole exon deletion of an out-of-frame exon which will ultimately result in restoration of the reading frame. Yet, the other criteria on the size of the skipped area (now deletion + skipped exon) and functional domains apply (see section titled “Notes on multi-exon Skipping”)</li>
        <li>It is theoretically possible to skip consecutive out-of-frame exons so as not to disrupt the reading frame, though designing an ASO that skips two consecutive exons may prove to be a challenge.</li>
        <li>It is theoretically possible to skip the first coding exon if a nearby in-frame start codon exists, and the first exon meets exon skipping criteria outlined in Table 5</li>
    </ul>

    <h5>Assessing Exon Size</h5>
    <p>
        Assessors should consider the size of the exon. As per ClinGen recommendations, an in-frame deletion in the size of 10% or more of the coding transcript is considered a strong criterion for loss of protein function (Walker et al., 2023). However, this is protein and exon-dependent (i.e., losing up to 30% of the dystrophin transcript can still result in a functional, truncated protein (Duan, 2016; Gao & McNally, 2015)). For this reason, skipping an exon that encodes for more than 10% of the protein is considered “unlikely eligible”, and obtaining functional evidence is the expected next step.
    </p>
    <p>
        The percentage of coding region can be calculated as follows:
    </p>
    <div style="background: #eee; padding: 10px; border-radius: 4px; font-family: monospace; margin: 1em 0;">
        Exon size as coding region in % = (exon length in bp/3) / length of the protein in aa *100
    </div>
    <p>
        Please see the example for the NM_003793.4(CTSF):c.264del (p.Cys89fs) [MIM 603539] variant in example video 10. The protein is 484 aa in length. Exon 2 is in-frame 0-0 and spans 99 nucleotides. This information can be obtained from Ensembl, UCSC genome browser, and/or UniProt (see video).
    </p>
    <p>
        Exon size in coding region % = (99/3) aa / 484 aa * 100 = 6.8 %<br>
        Alternatively, assessors can use the Genetic Therapy Generator Toolkit from the Dutch Center for RNA Therapeutics that automatically calculates the exon length for you. By submitting the variant in question, the GTGT will indicate how much protein coding region is left following skipping of the exon.
    </p>
    <p>
        There are instances where an exon coding for more than 10% of the protein should be classified as “not eligible”. If the exon encodes for more than 10% of the protein and codes for more than one non-repeat protein domain, this exon should be considered non-skippable and therefore not eligible for exon skipping ASOs. For more detail on assessing the role of functional domains, please see the section titled “Assessing the Role of Functional Domains”.
    </p>
    <p>
        In addition to the relative size (coding percentage), it is also important to consider an exon’s absolute size in nucleotides, as even large exons may constitute only a small portion of a very large protein. Large exons are generally more challenging therapeutic targets for ASO-mediated skipping because they often contain multiple splicing regulatory elements and cryptic splice sites, in theory resulting in the need for multiple ASOs to achieve effective exon skipping. An example is an ASO targetting exon 9 or the ATXN3 gene which ended up promoting a partial exon 9 skip via the activation of an alternative 5′ cryptic splice site (Evers et al., 2013). Instead, current successful exon skipping ASOs target moderately sized exons typically less than 300 nt, such as DMD exons 45, 51, and 53. A length less than 300 nt is believed to facilitate exon definition and recognition by the spliceosome, minimizing the likelihood of additional splicing regulatory elements and internal cryptic splice sites (Bolisetty & Beemon, 2012). Due to this, it is recommended that assessors be cautious when considering ASO-mediated skipping of exons of 300 nt or more.
    </p>
    <p>
        To further refine exon eligibility toward skipping, assessors can also use the recently developed VulExMap webtool, which classifies exons as either resilient or vulnerable to aberrant splicing as a result of mutations (Holm et al., 2024). Vulnerable exons, associated with less enforced inclusion and sensitivity to splicing perturbations, may be more amenable to ASO-mediated skipping, whereas resilient exons are robustly included in transcripts and may be more difficult to skip.
    </p>

    <h5>Searching for Cases of (Natural) Occurrences of Exon Skipping or In-frame Deletions</h5>
    <p>
        Assessors should determine whether exon skipping has already been observed naturally. This information is necessary to further determine whether exon skipping is a suitable option for a given variant. There are different ways in which exon skipping can occur naturally and resources such as ClinVar, gnomAD, DECIPHER, ExonSkipDB, and PubMed aid with their identification:
    </p>
    <p><strong>Canonical splice site or splice region variants</strong></p>
    <p>
        Assessors should search ClinVar, PubMed, gnomAD, or DECIPHER to collect data on variants that cause splice aberration leading to full exon skipping. This information should have been validated by sufficient functional data (e.g., RNAseq or qPCR). Note that variants affecting splicing motifs outside of the splice sites can also cause exon skipping and would also fall under this category.
    </p>
    <p>
        If full exon skipping has been observed and validated and is disease-causing (pathogenic variant), this exon (and therefore variants within that exon) is not eligible for exon skipping therapy. Please note that in silico predictions are not to be used as substitutes for RNAseq and qPCR data, and assessors must pay careful attention to what tools were used in the assessment of splicing outcomes of canonical splice site variants.
    </p>
    <p>
        If full-length exon skipping has been observed in individuals who do not show signs of the disease in question (individuals can of course have other diseases) the exon is eligible for exon skipping. If no evidence of pathogenic or benign exon skipping validated by qPCR or RNAseq exists, assessors should proceed with the analysis.
    </p>
    <p><strong>Full exon deletions or in-frame deletions within the exon</strong></p>
    <p>
        Assessors should search ClinVar, PubMed, gnomAD, or DECIPHER to collect data on in-frame deletions of full-exon deletions within the exon of interest. If in-frame deletions have been observed and validated to be disease-causing (pathogenic), this exon (and therefore variants within that exon) is not eligible. If this is the case, no further analysis is required. However, assessors should pay attention to whether the in-frame deletion creates a stop codon or new amino acid, as these exons could still be eligible for exon skipping via ASO (i.e., the pathogenicity is possibly a result of premature termination or the change of an amino acid leading to folding changes and not necessarily the deletion of amino acids itself). If this is the case, assessors should proceed with analysis.
    </p>
    <p>
        If full exon deletions and in some instances larger in-frame deletions are identified in the general population and classified as benign, this exon can be considered “eligible” for an exon skipping approach.
    </p>
    <p><strong>The assessment for the occurrence of natural exon skipping can easily be summarized as follows:</strong></p>
    <ul>
        <li>Eligible - Loss of the exon, either due to splice aberration or genomic deletions, has been classified as benign</li>
        <li>Not eligible - Loss of the exon or larger parts of the exon have been classified as pathogenic</li>
    </ul>

    <h5>Assessing the Role of Functional Domains</h5>
    <p>
        Assessors should also consider the functionality for which the exon codes. To begin, assessors should identify which domains are coded by the exon using tools such as the UCSC Genome Browser. Please cross-check domain annotations with UniProt or other sources to ensure the correct protein isoform is tagged. In this regard, we consider all types of domains as domains, including transmembrane and cytoplasmic domains as well as disordered regions. At this point, we do not have enough information and data to distill which domains can be removed without negative impact, yet, there are certain types of domains for which we can assume they are more important than others and of such high importance that skipping is not an option. Assessors can search for the role of the amino acids or domains present in an exon in the literature using PubMed, or utilize databases such as the Protein Database (PDB) or UniProt. Assessors should pay attention to the role of these specific domains (i.e., DNA binding domain, catalytic domains) or amino acids (i.e., specific amino acids known to play an important role in enzyme activity or protein structure). Sometimes, a web-search with “[protein name] protein structure” will also yield the desired results.
    </p>
    <p>
        Though it is difficult to define the importance of the domain or to predict the effects of exon skipping on protein function, assessors can consider the following exons as “not eligible” for exon skipping:
    </p>
    <ul>
        <li>The exon codes for important amino acid(s) with a known key functional role in the protein, or a functionally validated domain (i.e., involved in a catalytic domain, dimerization domain, inhibitory domain, etc.)</li>
        <li>The exon is a known mutational hotspot for pathogenic LoF missense variants (tools such as MetaDome, Essential3D, and Franklin can be used to help determine this). There is no strict definition of a mutational hotspot, but one can look at the distribution of the variants across the exons and specifically pay attention to the number of missense variants causing loss of protein function in a specific region/exon as this is highly indicative of an important functional domain. Further, if the ratio of LoF missense variants is higher than that of LoF truncating variants, this is a further indicator for an important functional domain (see also below “Special considerations for missense variants and small in-frame indels”). Web resources such as Franklin also have some calculations to indicate a mutational hotspot.</li>
        <li>The exon codes for the only functional domain in the protein (i.e., no possibility of residual function remaining)</li>
        <li>The exon codes for multiple functional domains and covers more than 10% of the coding region</li>
    </ul>
    <p>
        An exon that codes for a functional domain that does not meet the above criteria means it can be considered as “unlikely eligible”, indicating the need for functional studies. It is of note, however, that emerging research from newly developed ASOs has revealed exceptions to the above criteria. The ASO therapy zebronkysen, which promotes exon 8 skipping in the CLN3 (NM_001042432.2) gene to restore protein function as a potential treatment for Batten disease, has yielded promising preclinical and clinical results despite skipping an exon that encodes for a transmembrane domain (Cotman & Staropoli, 2012). This development reiterates the importance of reviewing literature on each exon of interest and the significance of functional studies following assessment. Nevertheless, it is still recommended that exons encoding for important functional domains not be considered high priority, and only functionally assessed if resource availability permits it.
    </p>
    <p>
        The above considerations apply whether the exon encodes part of a domain or the entire domain.
    </p>
    <p>
        Additionally, protein tandem repeat domains can also be considered as “unlikely eligible”, assuming the loss of a repeat (or part of a repeat domain) from the protein will still have a residual function (Duan, Goemans, Takeda, Mercuri, & Aartsma-Rus, 2021). Assessors can utilize the aforementioned strategies to assess tandem repeat domains. However, in cases where in-frame deletion of a repeat domain has been reported as pathogenic, or has been functionally proven to disrupt protein function, these exons would be considered “not eligible”. Most suitable are in-frame exons that contain a single, full repeat of the tandem repeats and no additional domains. For exon skipping of tandem repeat domains, we recommend that the protein consists of at least 5 tandem repeats. To identify exons that fulfill these criteria, assessors can utilize the tool TReXome DB
    </p>
    <p>
        Generally, we encourage contacting experts on a given gene and protein to discuss the role of the different domains and whether they consider exon skipping a possibility. Further, when looking into functional domains please check the exact role of a domain and the effect of losing the domain. For example, skipping an inhibitory domain could lead to a gain of function effect on protein level.
    </p>
    <p>
        If the exon codes for less than 10% of the protein, and no functional domain or important amino acids are coded for by the exon and the exon fulfills other criteria listed in Section A, then the exon can be considered as “likely eligible”. Exons are classified using the criteria outlined in Table 5.
    </p>

    <h5>Special considerations for missense variants and small in-frame indels</h5>
    <p>
        If a missense variant or small, in-frame indel leads to a LoF on the protein level, additional caution is necessary. These types of variants often indicate that the exon has a specific function even though no domain might yet be annotated. This could be that the respective amino acids are crucial for folding or the function of a certain domain has not yet been established. In these cases, we recommend paying special attention to other reported pathogenic variants within the exon. Should there be more pathogenic missense variants and in-frame indels than truncating variants, we consider this exon “unlikely eligible” for exon skipping.
    </p>
    <p>
        However, for missense variants within protein tandem repeat domains, we refer to the criteria outlined in the tandem repeat domain section above as this is an exception to this recommendation.
    </p>

    <h5>Notes on multi-exon skipping</h5>
    <p>
        Double- and multi-exon skipping have been suggested as approaches to bypass deleterious variants in out-of-frame exons when skipping adjacent exon(s) maintains the open reading frame. While this has been explored in preclinical studies, such as when targeting truncating variants in DMD exon 43 via 43-44 skipping (Aartsma-Rus et al., 2004), we currently do not recommend such an approach for clinical practice as it often requires the use of multiple ASOs.
    </p>
    <p>
        However, restoration of the reading frame in case of deletion of one or more out-of-frame exons can be achieved via skipping of a second out-of-frame exon that is adjacent to or downstream from a disease-causing exon deletion. If considering such an approach, criteria on exon size of the skipped area (including deleted + skipped exons), functional domains, and presence of missense variants/small in-frame indels still apply. In addition to this criteria, assessors should also consider utilizing tools such as GTGT, ExonViz, and the UCSC Genome Browser to identify instances where double- and multi-exon skipping produces an in-frame transcript. Assessors can also search for existing cases of double- or multi-exon deletions in the target region that are associated with mild to asymptomatic phenotypes as rationale for such an approach.
    </p>

    <h5>Important considerations for different inheritance patterns and pathomechanisms</h5>
    <p>
        The recommendations described above apply without restrictions to LoF variants in recessive disorders. For LoF variants in AD disorders, one might consider the development of allele-selectivespecific ASOs. This will ensure that exon skipping does not occur on the wildtype transcript and may result in greater amounts of functional gene product compared to non-allele selectivespecific ASOs. Though specific ASO design is beyond the scope of these guidelines, this may prove challenging and limit the types of ASOs to be designed.
    </p>
    <p>
        In the specific case of an AD disorder with an out-of-frame exon deletion, allele-selectivityspecificity is necessary otherwise exon skipping is not possible. Here, exon skipping of adjacent out-of-frame exons has the potential to restore the reading frame on the mutant allele, at the same time, exon skipping of that out-of-frame exon would destroy the reading frame of the wildtype allele.
    </p>
    <p>
        Further, for AD disorders associated with LoF variants, an upregulation of protein production from the wildtype allele should also be considered as an alternative option to exon skipping (see Section C).
    </p>
    <p>
        For GoF variants, additional considerations are possible. The guidelines for exon skipping can be applied as they are, with the exception that toxic GoF caused by the disruption of an inhibitory domain will not be rescued by skipping the exon containing the inhibitory domain (Mohassel et al., 2021). For GoF variants also an out-of-frame exon can be skipped, either the one containing the variant or any other one to downregulate transcript levels. In such a cases, haploinsufficiency of the gene should be taken into account which might require allele-selective approaches. Generally, for the downregulation of a transcript in the case of GoF and DN variants, Section B can be consulted. We do not recommend skipping out-of-frame exons for DN variants unless the WT allele remains intact (allele-selective skipping). For DN variants, it is important that the ration of WT to mutant allele needs to be shifted towards more WT allele.
    </p>

    <h5>ASO check</h5>
    <p>
        In addition to the recommended assessment strategies, assessors should review the literature or registries for exon skipping ASO strategies. This review can be performed either as the final step to validate your assessment strategy or earlier in the assessment process (see Step 0). Specifically, for exon skipping ASOs, it is crucial to evaluate whether an exon skipping approach has been implemented and validated at both the RNA and protein levels and shown to have the desired effect on the phenotype.
    </p>
    <p>
        In cases of conflicting evidence of ASO developments and the feasibility of an ASO therapy, consider the quality of the research, the types of experiments conducted, the nature of the results shared, and the publication date. The evaluation of literature on existing ASOs should be carried out at the assessor's discretion, with a critical and discerning approach.
    </p>
    <p>
        For exon skipping ASOs specifically, the approach will be applicable to all variants within a given exon. However, please pay attention to whether an existing ASO would bind to the variant site and could thus not be used for the individual under assessment if it is a different variant. Also check if the existing target site contains a SNP and wether the same SNP is present on the correct allele in the case under assessment. The exon (and therefore variant) is still considered “eligible” for canonical exon skipping according to the criteria outlined in Table 5 even if the ASO binds to a different variant site. This, however, means that a new ASO needs to be designed and developed for the variant under assessment.
    </p>
    <p>
        To help with identifying available ASOs for exon skipping, we recommend a search term in Pubmed like this, text in bold would need to be adjusted for the gene and exon being examined:
        <br><code>ABCA4 AND ((ASO) OR (AON) OR (antisense oligonucleotide)) AND (Exon 17)</code>
    </p>

    <h5>Table 5: Classification of variants for their eligibility towards exon skipping.</h5>
    <table class="classification-table">
        <thead>
            <tr>
                <th>Classification</th>
                <th>Criteria</th>
            </tr>
        </thead>
        <tbody>
            <tr class="bg-green">
                <td><strong>Eligible</strong></td>
                <td>Evidence that exon skipping does not impair protein function (benign canonical splice site variant leading to exon skipping, benign single exon deletion, naturally occurring transcript does not contain exon, or previously tested ASO with evidence of functional protein product).</td>
            </tr>
            <tr class="bg-blue">
                <td><strong>Likely eligible</strong></td>
                <td>
                    Exon is in-frame<br>
                    <strong>AND</strong><br>
                    Exon does NOT result in ≥10% loss of protein coding sequence if exon is skipped.<br>
                    <strong>AND</strong><br>
                    Exon does not create a stop codon OR new codon when skipped<br>
                    <strong>AND</strong><br>
                    Exon does not code for any functional domains<br>
                    <strong>AND</strong><br>
                    None of the exclusion criteria outlined in the “Not eligible” section are met.
                </td>
            </tr>
            <tr class="bg-orange">
                <td><strong>Unlikely eligible</strong></td>
                <td>
                    Exon is in-frame<br>
                    <strong>AND/OR</strong><br>
                    Exon does not create a stop codon when skipped BUT creates a new codon<br>
                    <strong>AND/OR</strong><br>
                    {<br>
                    Exon results in a loss of ≥10% of coding transcript if exon is skipped AND/OR exon codes for a single functional domain<br>
                    }<br>
                    <strong>AND</strong><br>
                    None of the exclusion criteria outlined in the “Not eligible” section are met
                </td>
            </tr>
            <tr class="bg-red">
                <td><strong>Not eligible</strong></td>
                <td>
                    Variant is in an out-of-frame exon OR<br>
                    Variant is in first or last coding exon OR<br>
                    Variant is in the ONLY coding exon OR<br>
                    Exon skipping results in a stop codon OR<br>
                    Exon skipping results in the loss of the ONLY functional domain in the protein OR<br>
                    Exon encodes for more than 10% of the proteins AND multiple non-repeat domains OR<br>
                    Exon codes for functionally proven important domains or amino acids (catalytic site, dimerization domain, inhibitory domains, etc.) OR<br>
                    Exon is a known mutational hotspot for (missense) loss-of-function variants OR<br>
                    Functional evidence of exon skipping shown to be pathogenic OR<br>
                    Functional evidence of in-frame deletions shown to be pathogenic OR<br>
                    Evidence that an ASO cannot be developed, shown by two independent investigations at the protein/functional level. Or one investigation with a convincing explanation why an ASO cannot be developed.
                </td>   
            </tr>
        </tbody>
    </table>
</div>

<div id="knockdown" class="tab-pane">
    <h4>Considerations for Transcript Knockdown</h4>
    <p>
        In this section, assessors will be guided on how to assess a variant for eligibility towards knockdown approaches. As described in the Background, knockdown ASOs/siRNAs can bind to the target transcript and downregulate (pre-)mRNA expression. Knockdown strategies can be utilized in cases where the pathomechanism is a result of overexpression, toxic GoF, or DN effects (Lauffer, van Roon-Mom, Aartsma-Rus, & N = 1 Collaborative, 2024).
    </p>

    <div class="toc-block" style="background: #f8f9fa; padding: 1em; border-left: 4px solid #005a9c; margin-bottom: 2em;">
        <strong>This section covers the following topics:</strong>
        <ul>
            <li>Considerations for pathomechanism</li>
            <li>Considerations for dosage sensitivity</li>
            <li>Important considerations for different inheritance patterns</li>
            <li>Strategies for searching the literature for gapmer ASOs and siRNAs</li>
        </ul>
    </div>

    <h5>Considerations for Pathomechanisms</h5>
    <p>
        Variants to be considered for knockdown approaches using gapmer ASOs or siRNA are toxic GoF and DN variants and copy number gains. Strategies on how to assess a variant mechanism are discussed in Step 2. Please ensure you have sufficient evidence of GoF or DN pathomechanism before proceeding further.
    </p>

    <h5>Considerations for Dosage Sensitivity</h5>
    <p>
        Before proceeding with the development of a knockdown ASO, it is crucial to consider dosage sensitivity. Ideally, a knockdown strategy would be employed when loss-of-function is not expected to cause disease. However, such cases are rare. More commonly, alterations in gene dosage are an underlying cause of disease and must be considered when developing knockdown ASOs. Here it becomes important to be aware of the different inheritance patterns that are implicated in a gene.
    </p>
    <p>
        If complete loss-of-function is not tolerated, but the loss of one gene copy is, knockdowns can still be considered. Resources such as PubMed, pLI scores, LOEUF scores (both available via gnomAD), DECIPHER dosage sensitivity track (using pHaplo, Collins et al., 2022), and ClinGen dosage sensitivity score can be used to determine this. We consider a curation from a reputable, independent source, such as the ClinGen consortium, the highest level of evidence for dosage sensitivity.
    </p>
    <p>
        An indication of whether the loss of one allele is tolerated can also be gained from population databases such as gnomAD; if there are carriers of LoF variants in the general population, it can be assumed that the loss of one allele is safe. The same considerations apply to carriers of homozygous LoF variants within a gene that implies that the complete loss of this gene is tolerated. Typically, if the loss of one gene copy is not a known cause of disease and has been observed in healthy cohorts, knockdown ASOs can be considered.
    </p>
    <p>
        Another factor to consider is haploinsufficiency. As described in Step 2, haploinsufficiency refers to a situation in which one healthy, wildtype allele does not generate sufficient protein product to preserve the physiological state (Deutschbauer et al., 2005). One can utilize resources such as gnomAD’s pLI and LOEUF scores, pHaplo score, and ClinGen’s dosage sensitivity score (preferred) to determine whether haploinsufficiency is a cause of disease. Typically, a pLI score of equal to or greater than 0.9, the top three deciles of LOEUF scores, a pHaplo score of 0.86 or greater, or a ClinGen haploinsufficiency score of 3 (sufficient evidence) indicates haploinsufficiency being a cause of disease (please note that this should be verified by functional evidence available in the literature). In cases where haploinsufficiency is a known cause of disease, a knockdown approach should only be considered if the disease caused by haploinsufficiency is less severe than the disorder associated with the GoF or DN variants. Though phenotype considerations are beyond the scope of this guideline, one can consider using OMIM, PubMed, and GeneReviews to assess the genotype-phenotype relation.
    </p>
    <p>
        One should also consider whether hypomorphic alleles are a cause of disease. Hypomorphic alleles are alleles which show partial loss-of-function (sometimes referred to as “leaky” alleles because there is some retention of protein function). In such cases, it is important to consider associated phenotypes, especially if partial loss-of-function is disease causing.
    </p>
    <p>
        For diseases that are tolerant to complete loss of function (knockout), targeting both alleles may be tolerated. The use of an allele-selectivespecific ASO targeting a SNP or the variant site should be considered if changes in gene dosage are a known cause of disease. Though this would greatly limit the ASO design, allele-selectivespecific ASOs in these scenarios would help ensure that at least 50% of the wildtype function remains (by targeting only the mutant allele).
    </p>
    <p>
        Another scenario in which allele-selectivespecific ASOs should be considered is if the mechanism is dominant-negative. DN variants impact the wildtype product and therefore may result in a functional product loss of greater than 50%. Hence, it is critical that an ASO is designed to specifically target the DN allele while keeping the wildtype product intact to recapitulate as much function as possible. For this reason, we recommend that, in the specific case of DN variants, that a variant can only be considered as likely eligible for a knockdown ASO only when an allele-selectivespecific approach can be taken.
    </p>

    <h5>Important considerations for different inheritance patterns</h5>
    <p>
        In the case of X-linked and Y-linked disorders, further considerations are necessary. For males, downregulation of a gene on the X or Y chromosome can lead to a complete loss of the gene function and before embarking on such an approach it should be known that this complete loss is tolerated. For X-linked disorders in females, it should be considered how X inactivation (discussed in Step 1) will affect gene dosage and whether a knockdown approach is safe. When considering a knockdown ASO, it is necessary to know whether the gene of interest escapes XCI or might reflect skewed XCI. Each of these situations will result in a different gene dosage, and therefore a different level of maximum expression to reduce.
    </p>
    <p>
        One X-linked case example is CLCN4-related condition, an XLD disease caused by a single pathogenic variant in the CLCN4 gene (Palmer, E. E. et al., 2018). It affects males and females differently depending on the pathomechanism (Palmer, Elizabeth Emma et al., 1993). LoF variants do no harm in females, but cause the disease in males. On the other hand, DN variants cause disease in females. This means that half the gene dose is sufficient, but anything less is not. An allele-selectivespecific knockdown ASO would therefore be required for DN variants in this gene for females only.
    </p>

    <h5>ASO Check</h5>
    <p>
        In addition to the recommended assessment strategies, assessors should review the literature for knockdown ASO/siRNA strategies. This review can be performed either as the final step to validate the assessment strategy or earlier in the assessment process (see Step 0). Specifically, for knockdown ASOs, it is crucial to evaluate whether a knockdown approach has been developed for other DN or GoF variants in the same gene, and that this approach has been validated at the RNA and protein level and shown to rescue the phenotype (pre-clinical work is sufficient). Take note if an allele-selectivespecific approach was used, as this may limit ASO design and affect outcomes. Further, an ASO might have been developed that is specific for a variant (other than the one under assessment) or a SNP. Information on phasing of an individual’s variant with that SNP will then have to be obtained. The gene is still “eligible” for knockdown according to the criteria outlined in Table 6 even if the already available ASO is not suitable for the individual under assessment. This means a new ASO would have to be designed and developed for that individual.
    </p>
    <p>
        In cases of conflicting evidence regarding ASO developments and the feasibility of an ASO therapy, consider the quality of the research, the types of experiments conducted, the nature of the results shared, and the publication date. The evaluation of literature on existing ASOs should be carried out at the assessor's discretion, with a critical and discerning approach.
    </p>
    <p>
        To help with identifying available ASOs/siRNA for knockdown approaches, we recommend a search term in Pubmed like this, text in bold would need to be adjusted for the gene and variant being examined:
    </p>
    <div style="background: #f1f1f1; padding: 10px; border-radius: 4px; font-family: monospace; margin-bottom: 1em;">
        SCN2A AND ((ASO) OR (AON) OR (antisense oligonucleotide) OR (AOs) OR (siRNA) OR (RNAi) OR (gapmer) or (knockdown))<br><br>
        SCN2A AND ((ASO) OR (AON) OR (antisense oligonucleotide) OR (AOs) OR (siRNA) OR (RNAi)) AND ((p.R853Q) OR (p.Arg853Gln) OR (c.2558G>A))
    </div>
    <p>
        Once all this information is collected, one can use Table 6 to classify the variant’s eligibility towards ASO knockdowns as “eligible”, “likely”, “unlikely”, or “not eligible”.
        Variants, where not enough evidence exists on the pathomechanism or the dosage sensitivity cannot be assessed until more evidence is collected.
    </p>

    <h5>Table 6: Classification of variants for their eligibility towards knockdown</h5>
    <table class="classification-table">
        <thead>
            <tr>
                <th>Classification</th>
                <th>Criteria</th>
            </tr>
        </thead>
        <tbody>
            <tr class="bg-green">
                <td><strong>Eligible</strong></td>
                <td>ASO/RNAi/siRNA has already been developed and shown to work with available functional evidence (i.e., evidence of knockdown rescuing function, pre-clinical data is sufficient)</td>
            </tr>
            <tr class="bg-blue">
                <td><strong>Likely eligible</strong></td>
                <td>
                    Variant is Gain-of-Function or Dominant-Negative (functionally proven)<br>
                    <strong>AND</strong><br>
                    {<br>
                    Gene is tolerant to the reduction of gene dosage (i.e., gene is NOT haploinsufficient)<br>
                    <strong>AND/OR</strong><br>
                    Individuals with heterozygous LoF variants are present in population databases/described in medical literature, such that high penetrance for probably severe disease phenotypes are unlikely<br>
                    }
                </td>
            </tr>
            <tr class="bg-orange">
                <td><strong>Unlikely eligible</strong></td>
                <td>
                    Variant is Gain-of-Function or Dominant-Negative (functionally proven)<br>
                    <strong>BUT</strong><br>
                    Heterozygous LoF/Haploinsufficiency/Hypomorphic variants has/have been associated with a disease<br>
                    <strong>OR</strong><br>
                    X-linked gene that undergoes X-chromosome inactivation
                </td>
            </tr>
            <tr class="bg-red">
                <td><strong>Not eligible</strong></td>
                <td>
                    Intolerant to reduction (i.e., gene dosage is tightly regulated in humans, and knockdown is expected to lead to serious phenotypic consequences)<br>
                    <strong>OR</strong><br>
                    Evidence that an ASO cannot be developed, shown by two independent investigations at the protein/functional level. Or one investigation with a convincing explanation why an ASO cannot be developed.
                </td>
            </tr>
        </tbody>
    </table>
</div>

<<div id="wt" class="tab-pane">
    <h4>Considerations for Upregulation from the Wildtype Allele</h4>
    <p>
        For disorders caused by haploinsufficiency, one functional wildtype gene copy remains. In these situations, one can use ASOs to upregulate the wildtype allele. This approach, also known as targeted augmentation of nuclear gene output (TANGO) (Lim et al., 2020), can include skipping of poison exons, downregulating naturally occurring antisense transcripts, and targeting UTR regulator elements such as upstream open reading frames, all of which can increase the gene product. For a more detailed explanation of these approaches, see the Background, Fig. 3, Fig. 11, and the following publications:
    </p>
    <ol>
        <li>Lim et al. (2020) <a href="https://doi.org/10.1038/s41467-020-17093-9" target="_blank">https://doi.org/10.1038/s41467-020-17093-9</a></li>
        <li>Mittal et al. (2022) <a href="https://doi.org/10.1016/j.medj.2022.08.006" target="_blank">https://doi.org/10.1016/j.medj.2022.08.006</a></li>
        <li>Felker et al. (2023) <a href="https://doi.org/10.1016/j.gim.2023.100884" target="_blank">https://doi.org/10.1016/j.gim.2023.100884</a></li>
    </ol>

    <p>
        For convenience, the data found in the supplementary tables and databases from the three aforementioned publications have also been compiled into one excel sheet (most recent excel sheet can be found on the N1C website). Assessors can use this sheet to look up poison exons, naturally occurring antisense transcripts, and upstream open reading frames found in each of these publications. It is highly encouraged that assessors do not rely on this one resource alone, but also utilize it in conjunction with the other strategies discussed below. Throughout this section, we will refer to tools and databases that can be utilized to better supplement and verify the findings from the three papers.
    </p>
    
    <div style="background-color: #f9f9f9; padding: 15px; border-left: 5px solid #ff9800; margin: 20px 0;">
        <p style="margin: 0;">
            Please note that, unlike the other ASO strategies, these guidelines will not provide details on how to classify variants as “likely”, “unlikely", or “not eligible”. The upregulation of wildtype transcripts through the aforementioned strategies heavily depends on the availability of functional evidence for key regulatory elements. These strategies are not comparable to one another, and different approaches must be employed depending on the availability of regulatory elements. Instead, readers can reference this section to learn about the different strategies and resources they can utilize for their own analyses. However, in the case where a wildtype upregulation approach has already been established for a given gene with sufficient functional evidence, we consider a variant as “eligible” for wildtype upregulation. The sole existence of a poison exon, NAT or other information in the papers and our combined table without further data does not qualify as sufficient evidence and does not support the classification “eligible” for wildtype upregulation.
        </p>
    </div>

    <div class="toc-block" style="background: #f8f9fa; padding: 1em; border-left: 4px solid #005a9c; margin-bottom: 2em;">
        <strong>This section covers:</strong>
        <ul>
            <li>Targeting naturally occurring antisense transcripts</li>
            <li>Targeting upstream open reading frames</li>
            <li>Targeting poison exons and non-productive alternative splicing events</li>
        </ul>
    </div>

    <h5>Naturally Occurring Antisense Transcripts</h5>
    <p>
        Through downregulating naturally occurring antisense transcripts via gapmer ASOs, it is - in some instances - possible to upregulate gene expression (Fig. 3, and Example 9 discussing upregulation of UBE3A (Dindot et al., 2023)). Assessors can utilize resources such as PubMed, UCSC Genome Browser, HUGO Gene Nomenclature Committee, and Ensembl to search for known antisense transcripts for a given gene. If searching through a genome browser like UCSC, naturally occurring antisense transcripts would be found to run in the reverse direction on the complementcompliment strand, and would be noted as a non-coding RNA. It is important to verify that this transcript is indeed a non-coding antisense transcript, and not another protein coding gene. Not all complimentary transcripts are NATs. To verify, we suggest utilizing databases such HUGO, paying attention to alias names or symbols that may indicate itsit’s function as a non-coding antisense transcript (i.e. “AS”, “antisense”, etc.). Please also note that antisense transcripts can be located in different positions in respect to the sense transcript, which will influence whether they are good targets. To learn more about antisense transcripts and their directions, we recommend to further study Werner et al. 2024 and Zinad et al 2017.
    </p>
    <p>
        Assessors should also consider the level at which these antisense transcripts are expressed in the tissue of interest. Databases such as GTEx and the Human Protein Atlas can assist with this. No general threshold of expression level can be provided, as this is gene dependent. However, assessors should pay attention to which tissues the transcript is enriched in, and how expression levels compare to the gene we wish to upregulate. Note that the existence of an antisense transcripts alone is not enough to consider an ASO development. Factors such as tissue expression, regulatory function, and orientations of the antisense transcripts (see above) must be well understood before considering this approach. For this reason, we highly encourage a literature search through PubMed to supplement the above assessment. Functional evidence must be available, describing the role of the transcript as a downregulator of the gene of interest. Search terms can include the name of the gene of interest, followed by the term “antisense transcript”, or the official name of the antisense transcript or host gene itself.
    </p>

    <h5>Upstream Open Reading Frames</h5>
    <p>
        Targetting uORFs for upregulation of the WT allele can be approached in two different ways. One option is that an ASO can be designed to target and block the regulatory elements in untranslated regions, including uORFs (Fig. 3). Resources such as PubMed can be used to check whether there is a uORF that can be targeted in the relevant transcript. Additionally, the supplemental table 2 in the paper by Mittal et al. (2022) indicates for which genes uORFs have been identified. Note that not all uORFs act through an inhibitory mechanism, and proper characterization of their mechanism is essential in determining their eligibility as ASO targets. To collect this information, we recommend a review of the literature and eventually, functional studies. The uORF should also be present in relevant and highly expressed transcripts in the target tissue.
    </p>
    <p>
        In addition to blocking uORFs (as described in Fig. 3), another approach is skipping an exon in the UTR that contains a uORF start site (Fig. 11). This strategy has been explored in a pre-print by Feng et al. (2025). The study identified several transcript isoforms of BDNF regulated by uORFs. Experimental validation confirmed uORF-mediated regulation in five BDNF transcripts. One transcript showed marked upregulation when an exon was removed from its 5′ UTR. Therefore, if considering this approach, one must assess the following:
    </p>
    <ul>
        <li>Expression levels of the alternate transcript regulated by the uORF (verify that it is expressed in the tissue of interest at sufficient levels, tools such as GTEx can be used for this)</li>
        <li>uORFs are experimentally validated</li>
        <li>A “skippable exon” in the context of pORF upregulation is identified, defined as 5’ UTR exons containing a minimum of one uORF start codon that does not contain the translation start site</li>
    </ul>

    <img src="{{ url_for('static', filename='fig11.png') }}" class="responsive-img" alt="Skipping of untranslated exons">
    <p class="img-caption">
        <strong>Figure 11: Skipping of untranslated exons to promote translation of pORF</strong><br>
        ASO can be used to skip untranslated exons containing start sites of the upstream open reading frame (uORF), promoting translation of the primary open reading frame (pORF), increasing wildtype gene product.
    </p>

    <h5>Poison Exons and Non-Productive Alternate Splicing Events</h5>
    <p>
        One can design ASOs which target non-productive alternate splicing events including poison exons (Fig. 3). These ASOs would utilize the splice-switching approaches described in the splice-correction and exon skipping sections to promote canonical splicing of the wildtype transcript. Resources such as PubMed, Ensembl, VastDB, and the UCSC Genome Browser can aid in determining alternate splicing events in the transcript of interest. Additionally, supplemental table 2 from Mittal et al. (2022), supplemental data 2 from Lim et al. (2020), and supplemental data 1 from Felker et al. (2023) list all identified poison exons in these papers. Databases such as GTEx can further assist by determining the expression level of alternate transcripts in the target tissue. As noted earlier, no specific threshold is available when reviewing transcript issue, as this will differ by gene and disease considerations. Note that some alternate splicing events are crucial for the production of important transcripts and isoforms. Generally, skipping a poison exon should lead to increased WT levels high enough to alter the disease state. The most prominent example of poison exon skipping in SCN1A for Dravet syndrome was possible because up to 80% of the transcripts contain this poison exon (Tang et al., 2025).
    </p>
    <p>
        For more resources and strategies on upregulating wildtype gene products, please reference the Useful Tools section.
    </p>

    <h5>Important considerations for different inheritance patterns and pathomechanisms</h5>
    <p>
        The considerations outlined here are mainly applicable to diseases associated with haploinsufficiency, which is most likely caused by LoF variants in AD disorders. In rare cases, one could also consider applying upregulation of wildtype allele strategies in X-linked disorders in females where there is sufficient evidence that upregulation from the second X chromosome is possible. Please always consider the challenges caused due to X inactivation.
    </p>

    <h5>ASO Check</h5>
    <p>
        Also for the upregulation of wildtype allele approaches, one can search whether specific strategies already exist, are under development, or are pursued in clinical trials. Multiple strategies can be utilized to upregulate gene products from the wildtype transcript, and it is therefore important to employ a variety of strategies/approaches in the search terms to ensure a comprehensive review of the field.
    </p>
    <p>
        While we do not classify the different upregulation approaches, we consider a variant as “eligible” for upregulation if an upregulation strategy has been developed and demonstrated to work with sufficient evidence pre-clincally.
    </p>
</div>

<script>
    function openTab(tabId) {
        // Hide all tab panes
        document.querySelectorAll('.tab-pane').forEach(el => el.classList.remove('active'));
        // Deactivate all buttons
        document.querySelectorAll('.tab-btn').forEach(el => el.classList.remove('active'));
        
        // Show specific tab and activate button
        document.getElementById(tabId).classList.add('active');
        // Find the button that called this function based on the onclick attribute matching logic is unnecessary hardcoded, 
        // simpler is to find the button by text or index, but let's iterate to find the matching onclick attribute for cleaner logic 
        // actually easier: simpler to just use event.currentTarget if passed, but to keep the signature simple inside the HTML:
        const buttons = document.querySelectorAll('.tab-btn');
        buttons.forEach(btn => {
            if(btn.getAttribute('onclick').includes(tabId)) {
                btn.classList.add('active');
            }
        });
        // Update URL hash for deep linking
        if (history.pushState) {
            history.pushState(null, null, '#' + tabId);
        }
    }
</script>
{% endblock %}
"""

    cite_html = """
{% extends "base.html" %}
{% block content %}
<h3>How to cite</h3>
<p>
    Information on how to cite will be available here as soon as AVEC is published.
</p>
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
rcnv_df: Optional[pd.DataFrame] = None

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

        # Load SSCVDB from Excel
        try:
            sscvdb_df = pd.read_excel(sscvdb_path)
            sscvdb_df.columns = sscvdb_df.columns.str.strip()
            if 'Variant ID' in sscvdb_df.columns:
                sscvdb_df['Variant ID'] = sscvdb_df['Variant ID'].astype(str).str.strip()
        except Exception as e:
            print(f"Warning: Could not load SSCVDB.xlsx: {e}")

        # Load N1C supplementary table (uORF / NAT / PE per gene)
        try:
            n1c_supp_df = pd.read_excel(n1c_supp_path)
            n1c_supp_df.columns = n1c_supp_df.columns.str.strip()
            for col in ['Gene', 'uORF', 'NAT', 'PE']:
                if col in n1c_supp_df.columns:
                    n1c_supp_df[col] = n1c_supp_df[col].astype(str).str.strip()

        except Exception as e:
            print(f"Warning: Could not load N1C_Variant_Supp_Table.xlsx: {e}")

        # Fetch and load N1C variants data
        response = requests.get(N1C_API_URL, timeout=30)
        response.raise_for_status() # Will raise an error if the request fails
        n1c_data = response.json()
        n1c_variants_df = pd.DataFrame(n1c_data)
        
        # Sanitize the N1C columns we will search on
        if 'Gene' in n1c_variants_df.columns:
            n1c_variants_df['Gene'] = n1c_variants_df['Gene'].astype(str).str.strip()
        if 'Coding DNA change (c.)' in n1c_variants_df.columns:
            n1c_variants_df['Coding DNA change (c.)'] = n1c_variants_df['Coding DNA change (c.)'].astype(str).str.strip()

        # Fetch and load N1C assessed variants (curated) data
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

        rcnv_path = os.path.join(DATA_DIR, 'rCNV.gene_scores.tsv')
        global rcnv_df
        if os.path.exists(rcnv_path):
            rcnv_df = pd.read_csv(rcnv_path, sep='\t')
            # Standardize gene column for matching
            if 'gene' in rcnv_df.columns:
                rcnv_df['gene'] = rcnv_df['gene'].astype(str).str.strip().str.upper()
        else:
            print(f"Warning: rCNV file not found at {rcnv_path}")

    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"An error occurred during database loading: {e}")
        exit(1)
# --- Add this method to EnsemblClient ---
def get_protein_metadata(self, protein_id):
    """Fetches metadata for a protein, including length."""
    if not protein_id: return None
    return self._get(f"/lookup/id/{protein_id}")

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
        domain_sources = {src.lower() for src in (
            'CDD','Pfam','SMART','PROSITE profiles','DisProt','PROSITE patterns','PRINTS','TIGRFAM','ProDom', 'mobidb-lite', 'RepeatsDB'
        )}
        repeat_disorder_sources = {
            'mobidb-lite'
        }
        allowed_sources = domain_sources | repeat_disorder_sources
        def _domain_url(feat: Dict[str, Any]) -> Optional[str]:
            src = str(feat.get('type', '')).strip().lower()
            identifier = str(feat.get('id') or feat.get('hit_id') or '').strip()
            interpro_id = str(feat.get('interpro') or '').strip()
            ident_upper, interpro_upper = identifier.upper(), interpro_id.upper()
            if src == 'pfam' and identifier:
                return f"https://www.ebi.ac.uk/interpro/entry/pfam/{ident_upper}/"
            if src == 'smart' and identifier:
                return f"https://www.ebi.ac.uk/interpro/entry/smart/{ident_upper}/"
            if src.startswith('prosite') and identifier:
                return f"https://prosite.expasy.org/{ident_upper}"
            if src == 'prints' and identifier:
                return f"https://www.ebi.ac.uk/interpro/entry/prints/{ident_upper}/"
            if src == 'tigrfam' and identifier:
                return f"https://www.ebi.ac.uk/interpro/entry/tigrfams/{ident_upper}/"
            if src == 'cdd' and identifier:
                return f"https://www.ncbi.nlm.nih.gov/Structure/cdd/cddsrv.cgi?uid={identifier}"
            if src == 'superfamily' and identifier:
                return f"https://www.ebi.ac.uk/interpro/entry/superfamily/{ident_upper}/"
            if src == 'disprot' and identifier:
                return f"https://disprot.org/{identifier}"
            if interpro_upper.startswith('IPR'):
                return f"https://www.ebi.ac.uk/interpro/entry/InterPro/{interpro_upper}/"
            return None
        def _display_name(feat: Dict[str, Any]) -> str:
            desc = feat.get('description') or feat.get('id') or feat.get('hit_id') or 'Domain'
            src = feat.get('type') or feat.get('source')
            return f"{desc} ({src})" if src else str(desc)
        def _is_domain_feature(feature):
            ftype = str(feature.get('type', '')).strip()
            desc = str(feature.get('description', '') or '').lower()
            return (
                ftype.lower() in allowed_sources
                or any(k in desc for k in ('disorder','repeat','coiled-coil','coiled coil','low complexity'))
            )
        preliminary_domains = []
        for f in all_features:
            if not _is_domain_feature(f): continue
            feat = dict(f)
            feat['source'] = feat.get('type')
            feat['url'] = _domain_url(feat)
            feat['display_name'] = _display_name(feat)
            preliminary_domains.append(feat)
        unique_domains = {}
        for feat in preliminary_domains:
            src_norm = str(feat.get('source') or feat.get('type') or '').lower()
            ftype_norm = str(feat.get('type') or '').lower()
            desc_norm = str(feat.get('description') or '').lower()
            id_norm = str(feat.get('id') or feat.get('hit_id') or '').lower()
            logic_norm = str(feat.get('logic_name') or '').lower()
            is_mobidb = any('mobidb' in val for val in (src_norm, ftype_norm, desc_norm, id_norm, logic_norm))
            is_disorder_like = any('disorder' in val for val in (ftype_norm, desc_norm, logic_norm, id_norm))
            is_mobidb_disorder = is_mobidb and is_disorder_like
            if is_mobidb_disorder:
                feat['source'] = 'MobiDBLite'
                feat['display_name'] = 'MobiDBLite disorder prediction'
                feat['id'] = feat.get('id') or feat.get('hit_id') or 'MobiDBLite'
            key = None
            if is_mobidb_disorder:
                key = ('mobidb-lite', 'disorder_prediction')
            else:
                key = feat.get('interpro') or (feat.get('type'), feat.get('start'), feat.get('end'), feat.get('description') or feat.get('id'))
            if is_mobidb_disorder and key in unique_domains:
                existing = unique_domains[key]
                merged = dict(existing)
                merged['start'] = min(existing.get('start', float('inf')), feat.get('start', float('inf')))
                merged['end'] = max(existing.get('end', float('-inf')), feat.get('end', float('-inf')))
                merged['description'] = existing.get('description') or feat.get('description') or "MobiDBLite disorder prediction"
                merged['display_name'] = existing.get('display_name') or "MobiDBLite disorder prediction"
                merged['id'] = existing.get('id') or feat.get('id') or feat.get('hit_id')
                unique_domains[key] = merged
                continue
            if key not in unique_domains: unique_domains[key] = feat
        return list(unique_domains.values())
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

    def lookup_symbol_expand(self, symbol):
        """Fetches expanded gene data for a given symbol."""
        data = self._get(f"/lookup/symbol/human/{symbol}", params={'expand': '1'})
        return data if isinstance(data, dict) else None
    
# --- Helper & Parsing Functions ---

def _evaluate_splice_variant_position(
    variant_hgvs: str, 
    vep_data: Dict[str, Any], 
    details: Dict[str, Any],
    exonic_pathogenic_user_input: Optional[str] = None,
    moa_user_input: Optional[str] = None
) -> Optional[Dict[str, Any]]:
    """
    Assesses a validated splice-altering variant based on its genomic position.
    This logic is shared by database-found variants and user-validated variants.
    """
    client = EnsemblClient()
    # Find the most relevant consequence to get gene and transcript
    all_consequences = vep_data.get('transcript_consequences', [])
    target_consequence = choose_best_consequence(all_consequences)
    if not target_consequence:
        return {"classification": "Unable to Assess", "reason": "Could not determine a target transcript from VEP response."}

    gene_symbol = target_consequence.get('gene_symbol')
    transcript_id = target_consequence.get('transcript_id')
    if not gene_symbol or not transcript_id:
        return {"classification": "Unable to Assess", "reason": "Could not determine gene symbol or transcript ID."}

    core_hgvs_match = re.search(r'(c\..*)', variant_hgvs, re.IGNORECASE)
    if not core_hgvs_match: 
        return None
    core_canonical_hgvs = core_hgvs_match.group(1).lower()

    result = {"details": details}
    consequence_terms = set(target_consequence.get('consequence_terms', []))
    
    is_intronic_by_consequence = 'intron_variant' in consequence_terms or 'splice_acceptor_variant' in consequence_terms or 'splice_donor_variant' in consequence_terms
    is_intronic_by_notation = '+' in core_canonical_hgvs or '-' in core_canonical_hgvs
    is_intronic = is_intronic_by_consequence or is_intronic_by_notation

    if is_intronic:
        dist_match = re.search(r'[+-](\d+)', core_canonical_hgvs)
        if not dist_match: 
            return None 
        dist = int(dist_match.group(1))
        result['details']['Distance to Canonical Splice Site'] = f"{dist} bp"
        
        if '+' in core_canonical_hgvs: # Downstream (e.g., c.123+1G>A)
            if dist <= 5:
                result.update({"classification": "Not Eligible", "reason": "Variant is a splice-altering variant located too close to the canonical splice site (<=+5bp). If there is functional evidence variant results in leaky exon skipping, consider exon inclusion."})
            elif 6 <= dist <= 50:
                result.update({"classification": "Unlikely Eligible", "reason": "Variant is a splice-altering variant located near the canonical splice site (+6-+50bp). AVEC only assesses dependent on location. If you have experimental validation that the branch point is not weakened, this should be considered 'likely eligible"})
            else:
                result.update({"classification": "Likely Eligible", "reason": "Variant is a splice-altering variant in a favorable deep-intronic position (>+50bp)."})
        
        elif '-' in core_canonical_hgvs: # Upstream (e.g., c.124-2A>G)
            if dist <= 5:
                result.update({"classification": "Not Eligible", "reason": "Variant is a splice-altering variant located too close to the canonical splice site (>=-5bp). If there is functional evidence variant results in leaky exon skipping, consider exon inclusion."})
            elif 6 <= dist <= 100:
                result.update({"classification": "Unlikely Eligible", "reason": "Variant is a splice-altering variant located near the canonical splice site (-6-(-100b)p). AVEC only assesses dependent on location. If you have experimental validation that the branch point is not weakened, this should be considered 'likely eligible'."})
            else:
                result.update({"classification": "Likely Eligible", "reason": "Variant is a splice-altering variant in a favorable deep-intronic position (<-100bp)."})
        
    elif any(c in consequence_terms for c in ['missense_variant', 'synonymous_variant', 'stop_gained', 'frameshift_variant', 'inframe_deletion', 'inframe_insertion']):
        # This is an exonic variant.
        if exonic_pathogenic_user_input is None:
            # Prompt the user for input on pathogenicity
            return {
                "classification": "Awaiting User Input",
                "reason": "This is an exonic splice-altering variant. To assess for splice correction, please confirm if the variant is pathogenic *without* its splicing effect.",
                "user_validation_prompt_exonic_pathogenicity": True,
                "moa_user_input": moa_user_input
            }
        elif exonic_pathogenic_user_input == 'yes':
            result.update({
                "classification": "Manual Assessment Required",
                "reason": "User confirmed the variant is pathogenic independent of its splicing effect. Manual lookup should be performed. Please refer to the about/methods splice correction page for alternative strategies."
            })
        else: # User answered 'no'
            # Assess based on distance to exon boundaries.
            gene_data = client.lookup_symbol_expand(gene_symbol)
            if not gene_data or 'Transcript' not in gene_data:
                 return {"classification": "Unable to Assess", "reason": f"Could not fetch expanded gene data for {gene_symbol}."}

            transcript_object = next((t for t in gene_data['Transcript'] if t.get('id') == transcript_id), None)
            if not transcript_object:
                return {"classification": "Unable to Assess", "reason": f"Could not find transcript {transcript_id} in expanded gene data."}
            
            all_exons = extract_exons_from_transcript(transcript_object)
            
            variant_pos = vep_data.get('start')

            if not variant_pos:
                return {"classification": "Unable to Assess", "reason": "Could not determine variant position for distance calculation."}

            all_distances = []
            for exon in all_exons:
                if exon.get('start') and exon.get('end'):
                    all_distances.append(abs(variant_pos - exon['start']))
                    all_distances.append(abs(variant_pos - exon['end']))
            
            if not all_distances:
                return {"classification": "Unable to Assess", "reason": "Could not calculate distances to any exon boundaries."}

            dist = min(all_distances)
            
            result['details']['Distance to Canonical Splice Site'] = f"{dist} bp"

            if dist <= 5:
                result.update({"classification": "Not Eligible", "reason": "User confirmed benign protein effect, but variant is within 5bp of an exon boundary. This may be eligible for exon inclusion if leaky exon skipping is observed."})
            elif 6 <= dist <= 15:
                result.update({"classification": "Unlikely Eligible", "reason": "User confirmed benign protein effect. Variant is between 6-15bp of an exon boundary."})
            else: # dist > 15
                result.update({"classification": "Likely Eligible", "reason": "User confirmed benign protein effect. Variant is >15bp from an exon boundary, making it a potential candidate for splice correction."})

    else:
        result.update({"classification": "Unable to assess", "reason": "Something went wrong."})

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

def _build_variant_link(variant: Dict[str, Any]) -> Optional[str]:
    """
    Builds a best-effort external link for a variant returned by Ensembl overlap endpoints.
    Prefers ClinVar links when available, otherwise falls back to dbSNP or Ensembl Variation explorer.
    """
    variant_id = str(variant.get('id') or variant.get('variation_name') or '').strip()
    if not variant_id:
        return None

    clinvar_id = variant.get('clinvar_variation_id') or variant.get('clinvar_id')
    if clinvar_id:
        cv_id = str(clinvar_id).strip()
        if cv_id:
            return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{cv_id}/"

    source = str(variant.get('source') or '').lower()
    if source == 'clinvar' and variant_id.isdigit():
        return f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{variant_id}/"

    if variant_id.lower().startswith('rs'):
        return f"https://www.ncbi.nlm.nih.gov/snp/{variant_id}"

    return f"https://www.ensembl.org/Homo_sapiens/Variation/Explore?v={variant_id}"

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

def _row_is_exon_skipping(row: pd.Series) -> bool:
    """Checks if a registry row describes an exon-skipping approach."""
    approach_cols = ['Approach', 'Therapeutic Modality', 'Therapeutic approach', 'Modality']
    for col in approach_cols:
        if col in row and isinstance(row.get(col), str):
            if re.search(r"exon\s*skip", row.get(col), re.IGNORECASE):
                return True
    # Fallback: search any string field for exon skip
    for col in row.index:
        val = row[col]
        if isinstance(val, str) and re.search(r"exon\s*skip", val, re.IGNORECASE):
            return True
    return False

def _row_is_knockdown(row: pd.Series) -> bool:
    """Checks if a registry row describes a knockdown approach."""
    approach_cols = ['Approach', 'Therapeutic Modality', 'Therapeutic approach', 'Modality']
    for col in approach_cols:
        if col in row and isinstance(row.get(col), str):
            if re.search(r"knock\s*down|silenc", row.get(col), re.IGNORECASE):
                return True
    for col in row.index:
        val = row[col]
        if isinstance(val, str) and re.search(r"knock\s*down|silenc", val, re.IGNORECASE):
            return True
    return False

def n1c_exon_skipping_exon_numbers_for_gene(gene_symbol: str) -> Tuple[set, List[str], Dict[int, List[str]]]:
    """
    Returns (set_of_exon_numbers, list_of_links, map_exon_to_links) for N1C registry rows that
    indicate exon skipping for the given gene. Filters to Approach: Exon Skipping when available.
    """
    exon_set: set = set()
    links: set = set()
    exon_link_map: Dict[int, set] = {}
    if n1c_variants_df is None or n1c_variants_df.empty or not gene_symbol:
        return exon_set, [], {}
    try:
        if 'Gene' in n1c_variants_df.columns:
            gene_matches = n1c_variants_df[n1c_variants_df['Gene'].astype(str).str.upper() == gene_symbol.upper()]
        else:
            gene_matches = n1c_variants_df
    except Exception:
        gene_matches = n1c_variants_df
    if gene_matches is None or len(gene_matches) == 0:
        return exon_set, [], {}

    candidate_cols = [
        'Target Exon', 'Targeted Exon', 'Exon', 'Exons', 'Exon(s)', 'Exon(s) Targeted',
        'Exon Target', 'Exon Targets', 'Intervention Description', 'Therapeutic Modality',
        'Approach', 'Description', 'Notes', 'Summary', 'Treatment Details'
    ]

    for _, row in gene_matches.iterrows():
        if not _row_is_exon_skipping(row):
            continue

        texts: List[str] = []
        for col in candidate_cols:
            if col in row and isinstance(row.get(col), str) and row.get(col).strip():
                texts.append(row.get(col))

        if not texts:
            try:
                texts = [val for val in row if isinstance(val, str)]
            except Exception:
                texts = []

        nums: List[int] = []
        for t in texts:
            nums.extend(_extract_exon_numbers_from_text(t))
        # Deduplicate in order
        seen = set()
        dedup_nums = []
        for n in nums:
            if n not in seen:
                dedup_nums.append(n)
                seen.add(n)

        if not dedup_nums:
            continue

        link = None
        try:
            nid = row.get('ID')
            if nid is not None and str(nid).strip() and str(nid).strip().lower() != 'nan':
                link = f"https://generegistry.n1collaborative.org/entry.html?id={str(nid).strip()}"
        except Exception:
            link = None

        for n in dedup_nums:
            exon_set.add(n)
            if link:
                links.add(link)
                exon_link_map.setdefault(n, set()).add(link)

    # Normalize map values to lists
    normalized_map: Dict[int, List[str]] = {k: sorted(v) for k, v in exon_link_map.items()}
    return exon_set, sorted(links), normalized_map

def n1c_exon_skipping_variant_exon_map(gene_symbol: str, all_exons: List[Dict[str, Any]], client: EnsemblClient) -> Dict[int, List[str]]:
    """
    Maps N1C registry exon-skipping variants (using RefSeq Transcript + c. change) to exon numbers,
    using the same HGVS parsing as user input. Returns {exon_number: [links]}.
    """
    exon_map: Dict[int, List[str]] = {}
    if n1c_variants_df is None or n1c_variants_df.empty or not gene_symbol or not all_exons:
        return exon_map

    try:
        if 'Gene' in n1c_variants_df.columns:
            gene_matches = n1c_variants_df[n1c_variants_df['Gene'].astype(str).str.upper() == gene_symbol.upper()]
        else:
            gene_matches = n1c_variants_df
    except Exception:
        gene_matches = n1c_variants_df

    if gene_matches is None or len(gene_matches) == 0:
        return exon_map

    for _, row in gene_matches.iterrows():
        if not _row_is_exon_skipping(row):
            continue

        refseq_val = None
        cdot_val = None
        for cand in ['RefSeq Transcript', 'Transcript', 'RefSeq', 'Transcript ID']:
            if cand in row and isinstance(row.get(cand), str) and row.get(cand).strip():
                refseq_val = row.get(cand).strip()
                break
        for cand in ['Coding DNA change (c.)', 'Variant (c.)', 'c_dot', 'c.']:
            if cand in row and isinstance(row.get(cand), str) and row.get(cand).strip():
                cdot_val = row.get(cand).strip()
                break

        if not refseq_val or not cdot_val:
            continue

        hgvs_input = f"{refseq_val}:{cdot_val}"
        try:
            hgvs_query, _ = parse_hgvs_query(hgvs_input)
            if not hgvs_query:
                continue
            vep_data = client.vep_hgvs(hgvs_query)
            if not vep_data or not isinstance(vep_data, list):
                continue
            vep_entry = vep_data[0]
            v_chrom, v_start, v_end = vep_entry.get('seq_region_name'), vep_entry.get('start'), vep_entry.get('end')
            if not all([v_chrom, v_start, v_end]):
                continue
            target_exon = next((ex for ex in all_exons if ex['seq_region_name'] == v_chrom and max(v_start, ex['start']) <= min(v_end, ex['end'])), None)
            if not target_exon:
                continue
            exon_num = target_exon.get('total_exon_number')
            if exon_num is None:
                continue
            link = None
            try:
                nid = row.get('ID')
                if nid is not None and str(nid).strip() and str(nid).strip().lower() != 'nan':
                    link = f"https://generegistry.n1collaborative.org/entry.html?id={str(nid).strip()}"
            except Exception:
                link = None
            if link:
                exon_map.setdefault(exon_num, []).append(link)
            else:
                exon_map.setdefault(exon_num, [])  # ensure key exists even without link
        except Exception:
            continue

    # Deduplicate links for each exon
    for ex_num, links in list(exon_map.items()):
        dedup = []
        seen = set()
        for l in links:
            if l not in seen:
                seen.add(l)
                dedup.append(l)
        exon_map[ex_num] = dedup
    return exon_map

def n1c_gene_knockdown_entry(gene_symbol: str) -> Optional[Dict[str, Any]]:
    """Returns a knockdown registry match at the gene level (Approach: Knockdown)."""
    if n1c_variants_df is None or n1c_variants_df.empty or not gene_symbol:
        return None
    try:
        if 'Gene' in n1c_variants_df.columns:
            gene_matches = n1c_variants_df[n1c_variants_df['Gene'].astype(str).str.upper() == gene_symbol.upper()]
        else:
            gene_matches = n1c_variants_df
    except Exception:
        gene_matches = n1c_variants_df
    if gene_matches is None or len(gene_matches) == 0:
        return None
    for _, row in gene_matches.iterrows():
        if not _row_is_knockdown(row):
            continue
        link = None
        try:
            nid = row.get('ID')
            if nid is not None and str(nid).strip() and str(nid).strip().lower() != 'nan':
                link = f"https://generegistry.n1collaborative.org/entry.html?id={str(nid).strip()}"
        except Exception:
            link = None
        modality = row.get('Therapeutic Modality') if isinstance(row.get('Therapeutic Modality'), str) else row.get('Approach')
        reason = (
            f"N1C registry lists knockdown project(s) for {gene_symbol}. "
            f"Because the variant is Autosomal Dominant with GoF mechanism, it is considered eligible for knockdown."
        )
        details = {}
        if modality:
            details["Registry Modality"] = str(modality)
        if link:
            details["N1C Registry Link"] = link
        return {
            "classification": "Eligible",
            "reason": reason,
            "link": link,
            "details": details
        }
    return None

def get_gene_characteristics(gene_symbol: str) -> Dict[str, Any]:
    """
    Retrieves MOI, Haploinsufficiency, Triplosensitivity, and MOA from loaded dataframes.
    Includes rCNV scores (pHaplo/pTriplo) from Collins et al. 2022.
    """
    characteristics = {
        "moi": [],
        "haploinsufficiency": {"text": "Unknown", "url": None},
        "triplosensitivity": {"text": "Unknown", "url": None}, # Added init
        "moa": [],
        "gene_url": None,
        "rcnv": {"pHaplo": "N/A", "pTriplo": "N/A", "url": "https://doi.org/10.1016/j.cell.2022.06.036"} # [NEW] Init rCNV
    }

    # Use the global rcnv_df loaded in load_databases()
    if rcnv_df is not None and gene_symbol:
        # Exact match (case-insensitive handled by load step)
        match = rcnv_df[rcnv_df['gene'] == gene_symbol.upper()]
        if not match.empty:
            ph = match.iloc[0].get('pHaplo')
            pt = match.iloc[0].get('pTriplo')
            # Format to 3 decimal places if numeric
            try:
                if pd.notna(ph): characteristics['rcnv']['pHaplo'] = f"{float(ph):.3f}"
                if pd.notna(pt): characteristics['rcnv']['pTriplo'] = f"{float(pt):.3f}"
            except (ValueError, TypeError):
                pass

    # --- ClinGen Lookup ---
    if clingen_df is not None and gene_symbol in clingen_df.index:
        gene_data_row = clingen_df.loc[gene_symbol]
        if isinstance(gene_data_row, pd.DataFrame):
            gene_data_row = gene_data_row.iloc[0]

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
            elif hap_score_str.startswith('40 -'): haplo_text = "Dosage sensitivity unlikely"
            else: haplo_text = "No evidence"
            
            characteristics['haploinsufficiency'] = {
                "text": haplo_text,
                "url": hap_url if pd.notna(hap_url) else None
            }

        # 3. Triplosensitivity with Link
        triplo_assertion = gene_data_row.get('dosage_triplosensitivity_assertion')
        triplo_url = gene_data_row.get('dosage_report')
        if pd.notna(triplo_assertion):
            triplo_score_str = str(triplo_assertion).strip()
            triplo_text = "Unknown"
            if triplo_score_str.startswith('3 -'): triplo_text = "Sufficient evidence"
            elif triplo_score_str.startswith('2 -'): triplo_text = "Emerging evidence"
            elif triplo_score_str.startswith('1 -'): triplo_text = "Little evidence"
            elif triplo_score_str.startswith('30 -'): triplo_text = "Gene associated with autosomal recessive phenotype"
            elif triplo_score_str.startswith('0 -'): triplo_text = "No evidence"
            elif triplo_score_str.startswith('40 -'): triplo_text = "Sensitivity unlikely"
            
            characteristics['triplosensitivity'] = {
                "text": triplo_text,
                "url": triplo_url if pd.notna(triplo_url) else None
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

def assess_knockdown(gene_characteristics: Dict[str, Any], overrides: Dict[str, bool] = {}) -> Dict[str, Any]:
    """Assesses eligibility for a knockdown strategy."""
    haplo_obj = gene_characteristics.get("haploinsufficiency", {"text": "Unknown"})
    haplo_status_text = haplo_obj.get("text", "Unknown")
    
    moi = gene_characteristics.get("moi", [])
    # --- NEW: Logic to track original vs overridden state ---
    check_name = "Gene is not haploinsufficient"
    original_passed = haplo_status_text in ["No evidence", "Little evidence", "Dosage sensitivity unlikely"]
    
    has_override = check_name in overrides
    final_passed = overrides.get(check_name, original_passed)

    checks = {
        check_name: {
            "passed": final_passed,
            "original_passed": original_passed,
            "overridden": has_override
        }
    }
    # Determine if haploinsufficiency status is unknown/missing for warning purposes
    haplo_unknown = (haplo_status_text is None) or (str(haplo_status_text).strip().lower() in ["", "unknown", "n/a"])
    
    # Classification logic now depends on the final 'is_not_haploinsufficient' value
    if not final_passed:
        return {
            "classification": "Unlikely Eligible", 
            "reason": "Gene is associated with haploinsufficiency or Autosomal Dominant LoF, therefore knockdown could lead to unintended consequences. If the haploinsufficiency phenotype is mild compared to the current phenotype, knockdown may still be considered.",
            "checks": checks
        }
    
    elif haplo_unknown and not has_override: # Only show 'Unable to Assess' if status is unknown AND not overridden
        reason = "Gene haploinsufficiency status is not available in datasources. Therefore knockdown could lead to unintended consequences and AVEC cannot determine amenability to knockdown."
        return {
            "classification": "Unable to Assess",
            "reason": reason,
            "checks": checks
        }
    else:
        # Base reason
        reason = "Gene is assessed as insensitive to haploinsufficiency by ClinGen or by manual assessment."

        return {
            "classification": "Likely Eligible",
            "reason": reason,
            "checks": checks
        }
    
def get_overlapping_genes(self, gene_id):
        """Fetches all genes that overlap with a given Ensembl Gene ID."""
        data = self._get(f"/overlap/id/{gene_id}", params={"feature": "gene"})
        return data if isinstance(data, list) else [] # This is a duplicate method

def assess_wt_upregulation(client: EnsemblClient, gene_id: str, gene_symbol: str) -> Dict[str, Any]:
    """
    Assesses for WT upregulation by checking for overlapping NATs and by
    searching for a conventional antisense gene name ([GENE]-AS1).
    Trusts the '-AS1' naming convention without a strict biotype check.
    """
    if not gene_id or not gene_symbol:
        return {"classification": "Unable to Assess", "reason": "Missing Gene ID or Symbol."}

    found_antisense_genes = {} 

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

        # Start details with curated statuses and cite the original guideline table
        details = dict(supp_details)
        details["Datasource: TableS2 of N1C VARIANT Guidelines"] = "https://www.sciencedirect.com/science/article/pii/S0002929725000643#app2"
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
                "classification": "Potential possibilities identified",
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
                "classification": "No potential possibilities identified",
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

def assess_splice_switching(variant_hgvs: str, vep_data: Dict[str, Any], gene_symbol: str, exonic_pathogenic_user_input: Optional[str] = None) -> Optional[Dict[str, Any]]:
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
        # Not found in SpliceVarDB for this gene -> check SSCVDB fallback
        if sscvdb_df is not None and not sscvdb_df.empty and vep_data:
            variant_key = _format_sscvdb_variant_id_from_vep(vep_data)
            if variant_key and 'Variant ID' in sscvdb_df.columns:
                if not sscvdb_df[sscvdb_df['Variant ID'].str.strip().str.lower() == variant_key.strip().lower()].empty:
                    details = {
                        "Source Database": "SSCVDB",
                        "Evidence": "Splice-altering reported in SSCVDB",
                        "SSCVDB Gene Page": f"https://sscvdb.io/gene/{gene_symbol}"
                    }
                    return _evaluate_splice_variant_position(variant_hgvs, vep_data, details, exonic_pathogenic_user_input=exonic_pathogenic_user_input)
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

            return _evaluate_splice_variant_position(variant_hgvs, vep_data, details, exonic_pathogenic_user_input=exonic_pathogenic_user_input)

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
                return _evaluate_splice_variant_position(variant_hgvs, vep_data, details, exonic_pathogenic_user_input=exonic_pathogenic_user_input)

    # --- Return prompt object if still not found ---
    return {
        "classification": "Not in Database",
        "reason": "This variant was not found in SpliceVarDB. If there is experimental validation (by qPCR or RNA-seq) that this variant is splice-altering, please confirm below.",
        "user_validation_prompt": True
    }

def assess_consecutive_exons(client, transcript, all_exons, target_exon, vep_entry, refseq_id_for_viewer, overrides={}):
    coding_exons = sorted([e for e in all_exons if e['cds_length'] > 0], key=lambda x: x['coding_exon_number'])
    
    target_idx = next((i for i, ex in enumerate(coding_exons) if ex['exon_id'] == target_exon['exon_id']), -1)
    if target_idx == -1: return None

    # Helper to assess a specific pair of exons
    def evaluate_pair(exon_a, exon_b):
        gene_id = transcript.get('Parent')
        transcript_id = transcript.get('id')
        protein_id = transcript.get("Translation", {}).get("id")
        protein_meta = client._get(f"/lookup/id/{protein_id}") if protein_id else None
        
        total_cds_len = sum(e['cds_length'] for e in coding_exons)
        total_protein_len = protein_meta['length'] if (protein_meta and 'length' in protein_meta) else (total_cds_len / 3 if total_cds_len > 0 else 0)
        
        combined_cds_len = exon_a['cds_length'] + exon_b['cds_length']
        combined_aa_len = combined_cds_len / 3
        combined_len_nt = (exon_a['end'] - exon_a['start'] + 1) + (exon_b['end'] - exon_b['start'] + 1)
        
        # Stop codon check
        cds_seq = client.get_cds_sequence(transcript_id)
        orig_cond2_no_stop = False
        orig_cond10_no_new_codon = True
        if cds_seq:
            try:
                current_pos = 0
                cds_map = {}
                for ex in coding_exons:
                    cds_map[ex['coding_exon_number']] = cds_seq[current_pos : current_pos + ex['cds_length']]
                    current_pos += ex['cds_length']
                skipped_cds = "".join(cds_map[i] for i in sorted(cds_map.keys()) if i != exon_a['coding_exon_number'] and i != exon_b['coding_exon_number'])
                if skipped_cds:
                    prot = str(Seq(skipped_cds).translate(to_stop=False))
                    orig_cond2_no_stop = "*" not in prot[:-1]
                
                # Check for novel codon creation at junction
                idx_a = next((i for i, ex in enumerate(coding_exons) if ex['coding_exon_number'] == exon_a['coding_exon_number']), -1)
                
                if idx_a > 0 and idx_a + 1 < len(coding_exons) - 1:
                    len_prev = sum(e['cds_length'] for i, e in enumerate(coding_exons) if i < idx_a)
                    phase = len_prev % 3
                    
                    if phase != 0:
                        prev_seq = cds_map[coding_exons[idx_a-1]['coding_exon_number']]
                        next_seq = cds_map[coding_exons[idx_a+2]['coding_exon_number']]
                        target_seq = cds_map[exon_a['coding_exon_number']]
                        needed = 3 - phase
                        if len(prev_seq) >= phase and len(next_seq) >= needed and len(target_seq) >= needed:
                            new_codon = prev_seq[-phase:] + next_seq[:needed]
                            orig_codon = prev_seq[-phase:] + target_seq[:needed]
                            if str(Seq(new_codon).translate()) != str(Seq(orig_codon).translate()):
                                orig_cond10_no_new_codon = False
            except Exception: orig_cond2_no_stop = False
            
        orig_cond3_not_terminal = (exon_a['coding_exon_number'] != 1 and exon_b['coding_exon_number'] != len(coding_exons))
        orig_cond4_small = (Decimal(str(combined_aa_len)) / Decimal(str(total_protein_len))) < Decimal('0.1') if total_protein_len > 0 else False
        
        domains = client.get_domains(protein_id) if protein_id else []
        overlapping_domain_names, overlapping_domain_details = [], []
        if domains:
            cds_pos_start = sum(e['cds_length'] for e in coding_exons if e['coding_exon_number'] < exon_a['coding_exon_number'])
            exon_aa_start = (cds_pos_start // 3) + 1
            exon_aa_end = ((cds_pos_start + combined_cds_len - 1) // 3) + 1
            for d in domains:
                if d.get('start', 0) <= exon_aa_end and d.get('end', 0) >= exon_aa_start:
                    label = d.get('display_name') or d.get('description', d.get('id', 'Unknown Domain'))
                    overlapping_domain_names.append(label)
                    overlapping_domain_details.append({
                        "label": label, "name": d.get('description') or d.get('id') or "Domain",
                        "source": d.get('source') or d.get('type'), "id": d.get('id'), "url": d.get('url'),
                        "start": d.get('start'), "end": d.get('end')
                    })
        domain_count = len(overlapping_domain_names)
        orig_cond5_no_domain = domain_count == 0
        
        # Variants
        def get_counts(chrom, start, end):
            vars_in_reg = client.overlap_region_variation(chrom, start, end)
            c = {'missense': 0, 'inframe_del': 0, 'splice': 0, 'nonsense': 0, 'frameshift': 0, 'benign_splice':0}
            links = {"pathogenic": [], "benign": []}
            seen = {"pathogenic": set(), "benign": set()}
            for v in vars_in_reg:
                clclass = classify_variant_clinsig(v.get('clinical_significance'))
                conseq = (v.get("consequence_type") or "").lower()
                vid = str(v.get('id') or v.get('variation_name') or '').strip()
                vsrc = (v.get('source') or '').strip() or "Ensembl"
                vlink = _build_variant_link(v) if vid else None
                if clclass == "pathogenic":
                    if "missense" in conseq: c['missense'] += 1
                    elif "inframe_deletion" in conseq: c['inframe_del'] += 1
                    elif "splice_donor" in conseq or "splice_acceptor" in conseq: c['splice'] += 1
                    elif "stop_gained" in conseq: c['nonsense'] += 1
                    elif "frameshift" in conseq: c['frameshift'] += 1
                    if ("splice_donor" in conseq or "splice_acceptor" in conseq) and vid and vid not in seen["pathogenic"]:
                        links["pathogenic"].append({"id": vid, "source": vsrc, "url": vlink})
                        seen["pathogenic"].add(vid)
                elif clclass == 'benign' and ("splice_donor" in conseq or "splice_acceptor" in conseq):
                    c['benign_splice'] += 1
                    if vid and vid not in seen["benign"]:
                        links["benign"].append({"id": vid, "source": vsrc, "url": vlink})
                        seen["benign"].add(vid)
            return c, links

        c1, l1 = get_counts(exon_a['seq_region_name'], exon_a['start'], exon_a['end'])
        c2, l2 = get_counts(exon_b['seq_region_name'], exon_b['start'], exon_b['end'])
        counts = {k: c1.get(k, 0) + c2.get(k, 0) for k in c1}
        splice_variant_links = {"pathogenic": l1["pathogenic"] + l2["pathogenic"], "benign": l1["benign"] + l2["benign"]}
        
        orig_cond6_missense = counts['missense'] < 3 + counts['frameshift'] + counts['nonsense']
        orig_cond7_splice = counts['splice'] == 0
        orig_cond8_no_inframe_del = counts['inframe_del'] == 0
        
        # Apply Overrides
        cond2_no_stop = overrides.get('No New Stop Codon', orig_cond2_no_stop)
        cond3_not_terminal = overrides.get('Not First/Last Exon', orig_cond3_not_terminal)
        cond4_small = overrides.get('Is <10% of Protein', orig_cond4_small)
        cond5_no_domain = overrides.get('No Domain Overlap', orig_cond5_no_domain)
        cond6_missense = overrides.get('Low Missense Count', orig_cond6_missense)
        cond7_splice = overrides.get('No Pathogenic Splice Variants', orig_cond7_splice)
        cond8_no_inframe_del = overrides.get('No Pathogenic In-Frame Deletions', orig_cond8_no_inframe_del)
        cond10_no_new_codon = overrides.get('No New Codon Formation', orig_cond10_no_new_codon)

        classification, reason = "Undetermined", ""
        if not cond3_not_terminal: classification, reason = "Not Eligible", "One of the exons is the first or last coding exon."
        elif not cond2_no_stop: classification, reason = "Not Eligible", "Skipping these exons is predicted to create a premature stop codon."
        elif not cond7_splice: classification, reason = "Not Eligible", f"Exons contain {counts['splice']} pathogenic splice variant(s)."
        elif not cond4_small and domain_count > 1: classification, reason = "Not Eligible", f"Exons overlap with {domain_count} protein domains and constitute >10% of the coding region."
        elif not cond8_no_inframe_del: classification, reason = "Not Eligible", f"Exons contain {counts['inframe_del']} pathogenic in-frame deletion(s)."
        elif not cond10_no_new_codon: classification, reason = "Unlikely Eligible", "Skipping these exons creates a novel codon at the junction that codes for a different amino acid."
        elif not cond6_missense: classification, reason = "Unlikely Eligible", f"Exons are a mutational hotspot with {counts['missense']} pathogenic missense variants."
        elif not cond4_small: classification, reason = "Unlikely Eligible", "Exons constitute >=10% of the protein."
        elif not cond5_no_domain: classification, reason = "Unlikely Eligible", f"Exons overlap with {domain_count} protein domain(s)."
        else: classification, reason = "Likely Eligible", "Exons meet the primary criteria for multi-exon skipping."
        
        note = f"Assessment based on skipping two consecutive exons: Exon {exon_a['total_exon_number']} and Exon {exon_b['total_exon_number']}."
        if combined_len_nt > 300:
            note += " Combined exon length exceeds 300 nucleotides. Proceed with caution."
        note += " Please note that some exceptions for exon skipping assessment apply. Please refer to the methods section for more details."
            
        visualization_data = None
        try:
            v_chrom, v_start, v_end = vep_entry.get('seq_region_name'), vep_entry.get('start'), vep_entry.get('end')
            padding = 1000
            visualization_data = {
                "locus": f"{v_chrom}:{max(1, v_start - padding)}-{v_end + padding}",
                "variantTrack": {"name": "Variant", "features": [{"chr": v_chrom, "start": v_start - 1, "end": v_end, "name": vep_entry.get('id', 'Variant')}]},
                "domainTrack": None
            }
        except: pass

        return {
            "classification": classification, "reason": reason, "note": note, "exon_length": combined_len_nt,
            "frac_cds": f"{((Decimal(str(combined_aa_len)) / Decimal(str(total_protein_len))) * 100):.2f}%" if total_protein_len > 0 else "N/A",
            "pathogenic_variant_counts": counts, "domain_count": domain_count, "domain_names": overlapping_domain_names,
            "domain_details": overlapping_domain_details, "splice_variant_links": splice_variant_links,
            "coding_exon_number": f"{exon_a['coding_exon_number']}+{exon_b['coding_exon_number']}",
            "total_exon_number": f"{exon_a['total_exon_number']}+{exon_b['total_exon_number']}",
            "gene_id": gene_id, "transcript_id": transcript_id,
            "clinvar_url": f"https://www.ncbi.nlm.nih.gov/clinvar/?term=GRCh38%3A{exon_a['seq_region_name']}%3A{min(exon_a['start'], exon_b['start'])}-{max(exon_a['end'], exon_b['end'])}",
            "checks": {
                "Is In-Frame": {"passed": True, "original_passed": True, "overridden": 'Is In-Frame' in overrides},
                "No New Stop Codon": {"passed": cond2_no_stop, "original_passed": orig_cond2_no_stop, "overridden": 'No New Stop Codon' in overrides},
                "Not First/Last Exon": {"passed": cond3_not_terminal, "original_passed": orig_cond3_not_terminal, "overridden": 'Not First/Last Exon' in overrides},
                "No Pathogenic Splice Variants": {"passed": cond7_splice, "original_passed": orig_cond7_splice, "overridden": 'No Pathogenic Splice Variants' in overrides},
                "No Pathogenic In-Frame Deletions": {"passed": cond8_no_inframe_del, "original_passed": orig_cond8_no_inframe_del, "overridden": 'No Pathogenic In-Frame Deletions' in overrides},
                "No Domain Overlap": {"passed": cond5_no_domain, "original_passed": orig_cond5_no_domain, "overridden": 'No Domain Overlap' in overrides},
                "No New Codon Formation": {"passed": cond10_no_new_codon, "original_passed": orig_cond10_no_new_codon, "overridden": 'No New Codon Formation' in overrides},
                "Low Missense Count": {"passed": cond6_missense, "original_passed": orig_cond6_missense, "overridden": 'Low Missense Count' in overrides},
                "Is <10% of Protein": {"passed": cond4_small, "original_passed": orig_cond4_small, "overridden": 'Is <10% of Protein' in overrides}
            },
            "visualization": visualization_data
        }

    candidates = []
    if target_idx > 0:
        prev = coding_exons[target_idx - 1]
        if (prev['cds_length'] + target_exon['cds_length']) % 3 == 0:
            candidates.append((prev, target_exon))
    if target_idx < len(coding_exons) - 1:
        nxt = coding_exons[target_idx + 1]
        if (target_exon['cds_length'] + nxt['cds_length']) % 3 == 0:
            candidates.append((target_exon, nxt))
            
    if not candidates: return None
    
    results = []
    for ex_a, ex_b in candidates:
        results.append(evaluate_pair(ex_a, ex_b))
        
    def rank(res):
        c = res['classification'].lower()
        if 'not eligible' in c: return 0
        if 'unlikely' in c: return 1
        if 'likely' in c: return 2
        if 'eligible' in c: return 3
        return -1
        
    results.sort(key=rank, reverse=True)
    best_res = results[0]
    
    if len(results) > 1:
        best_rank = rank(best_res)
        other_res = results[1]
        other_rank = rank(other_res)
        if best_rank > 0 and best_rank == other_rank:
            note = best_res.get('note', '')
            note += f" Skipping Exon {other_res['total_exon_number']} is also {other_res['classification']}."
            best_res['note'] = note.strip()
            
    return best_res

def assess_single_exon(client, original_query, transcript, all_exons, target_exon, vep_entry: Dict[str, Any], overrides: Dict[str, bool] = {}, refseq_id_for_viewer: Optional[str] = None):
    # --- Step 1: Data Gathering and Calculations ---
    # Set precision for Decimal calculations
    getcontext().prec = 10

    gene_id = transcript.get('Parent') 
    transcript_id = transcript.get('id')
    protein_id = transcript.get("Translation", {}).get("id")
    protein_meta = client._get(f"/lookup/id/{protein_id}") if protein_id else None
    
    cds_seq = client.get_cds_sequence(transcript_id)
    
    coding_exons = [e for e in all_exons if e['cds_length'] > 0]
    total_coding_exons = len(coding_exons)
    total_cds_len = sum(e['cds_length'] for e in coding_exons)
    if protein_meta and 'length' in protein_meta:
        total_protein_len = protein_meta['length']
    else:
        total_protein_len = total_cds_len / 3 if total_cds_len > 0 else 0

    if not target_exon.get('coding_exon_number'):
        return {"classification": "Unable to Assess", "reason": f"The variant maps to exon {target_exon['total_exon_number']}, which is non-coding."}
    
    chrom, start, end = target_exon['seq_region_name'], target_exon['start'], target_exon['end']
    variants_in_region = client.overlap_region_variation(chrom, start, end)
    
    clinvar_url = f"https://www.ncbi.nlm.nih.gov/clinvar/?term=GRCh38%3A{chrom}%3A{start}-{end}"

    counts = {'missense': 0, 'inframe_del': 0, 'splice': 0, 'nonsense': 0, 'frameshift': 0, 'benign_splice':0}
    splice_variant_links = {"pathogenic": [], "benign": []}
    seen_splice_variants = {"pathogenic": set(), "benign": set()}
    for v in variants_in_region:
        clclass = classify_variant_clinsig(v.get('clinical_significance'))
        conseq = (v.get("consequence_type") or "").lower()
        variant_id = str(v.get('id') or v.get('variation_name') or '').strip()
        variant_source = (v.get('source') or '').strip() or "Ensembl"
        variant_link = _build_variant_link(v) if variant_id else None
        if clclass == "pathogenic":
            if "missense" in conseq: counts['missense'] += 1
            elif "inframe_deletion" in conseq: counts['inframe_del'] += 1
            elif "splice_donor" in conseq or "splice_acceptor" in conseq: counts['splice'] += 1
            elif "stop_gained" in conseq: counts['nonsense'] += 1
            elif "frameshift" in conseq: counts['frameshift'] += 1
            if ("splice_donor" in conseq or "splice_acceptor" in conseq) and variant_id and variant_id not in seen_splice_variants["pathogenic"]:
                splice_variant_links["pathogenic"].append({
                    "id": variant_id,
                    "source": variant_source,
                    "url": variant_link
                })
                seen_splice_variants["pathogenic"].add(variant_id)
        elif clclass == 'benign' and ("splice_donor" in conseq or "splice_acceptor" in conseq):
            counts['benign_splice'] += 1
            if variant_id and variant_id not in seen_splice_variants["benign"]:
                splice_variant_links["benign"].append({
                    "id": variant_id,
                    "source": variant_source,
                    "url": variant_link
                })
                seen_splice_variants["benign"].add(variant_id)
    
    exon_cds_len = target_exon['cds_length']
    exon_aa_len = exon_cds_len / 3
    coding_exon_number = target_exon['coding_exon_number']
    
    # --- Step 1.5: Apply Overrides ---
    overridden_checks = {}
    has_override = False

    # --- Step 2: Condition Checks ---
    # These are the automatically calculated states
    original_conds = {}
    original_conds['cond1_inframe'] = (exon_cds_len % 3 == 0)
    original_conds['cond2_no_stop'] = False
    original_conds['cond10_no_new_codon'] = True
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
                original_conds['cond2_no_stop'] = "*" not in prot[:-1]
            
            # Check for novel codon creation at junction
            target_idx = next((i for i, ex in enumerate(sorted_coding_exons) if ex['coding_exon_number'] == coding_exon_number), -1)
            if original_conds['cond1_inframe'] and target_idx > 0 and target_idx < len(sorted_coding_exons) - 1:
                len_prev = sum(e['cds_length'] for i, e in enumerate(sorted_coding_exons) if i < target_idx)
                phase = len_prev % 3
                if phase != 0:
                    prev_seq = cds_map[sorted_coding_exons[target_idx-1]['coding_exon_number']]
                    next_seq = cds_map[sorted_coding_exons[target_idx+1]['coding_exon_number']]
                    target_seq = cds_map[coding_exon_number]
                    needed = 3 - phase
                    if len(prev_seq) >= phase and len(next_seq) >= needed and len(target_seq) >= needed:
                        new_codon = prev_seq[-phase:] + next_seq[:needed]
                        orig_codon = prev_seq[-phase:] + target_seq[:needed]
                        if str(Seq(new_codon).translate()) != str(Seq(orig_codon).translate()):
                            original_conds['cond10_no_new_codon'] = False
        except Exception: 
            original_conds['cond2_no_stop'] = False

    original_conds['cond3_not_terminal'] = (coding_exon_number is not None and coding_exon_number not in (1, total_coding_exons))
    original_conds['cond4_small'] = (Decimal(str(exon_aa_len)) / Decimal(str(total_protein_len))) < Decimal('0.1') if total_protein_len > 0 else False
    domains = client.get_domains(protein_id) if protein_id else []
    overlapping_domain_names, overlapping_domain_details = [], []
    if domains:
        cds_pos_start = sum(e['cds_length'] for e in sorted(coding_exons, key=lambda x: x['coding_exon_number']) if e['coding_exon_number'] < coding_exon_number)
        exon_aa_start, exon_aa_end = (cds_pos_start // 3) + 1, ((cds_pos_start + exon_cds_len -1) // 3) + 1
        for d in domains:
            if d.get('start', 0) <= exon_aa_end and d.get('end', 0) >= exon_aa_start:
                label = d.get('display_name') or d.get('description', d.get('id', 'Unknown Domain'))
                overlapping_domain_names.append(label)
                overlapping_domain_details.append({
                    "label": label,
                    "name": d.get('description') or d.get('id') or d.get('hit_id') or "Domain",
                    "source": d.get('source') or d.get('type'),
                    "id": d.get('id') or d.get('hit_id'),
                    "url": d.get('url'),
                    "start": d.get('start'),
                    "end": d.get('end')
                })

    domain_count = len(overlapping_domain_names)
    original_conds['cond5_no_domain'] = domain_count == 0
    original_conds['cond6_missense'] = counts['missense'] < 3 + counts['frameshift'] + counts['nonsense']
    original_conds['cond7_splice'] = counts['splice'] == 0
    original_conds['cond8_no_inframe_del'] = counts['inframe_del'] == 0
    original_conds['cond9_benign_splice'] = counts['benign_splice'] > 0

    # --- Step 2.5: Apply Overrides to get FINAL states ---
    # Start with original values, then overwrite if an override exists
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
    cond4_small = (Decimal(str(exon_aa_len)) / Decimal(str(total_protein_len))) < Decimal('0.1') if total_protein_len > 0 else False    
    domains = client.get_domains(protein_id) if protein_id else []
    overlapping_domain_names, overlapping_domain_details = [], []
    if domains:
        cds_pos_start = sum(e['cds_length'] for e in sorted(coding_exons, key=lambda x: x['coding_exon_number']) if e['coding_exon_number'] < coding_exon_number)
        exon_aa_start, exon_aa_end = (cds_pos_start // 3) + 1, ((cds_pos_start + exon_cds_len -1) // 3) + 1
        for d in domains:
            if d.get('start', 0) <= exon_aa_end and d.get('end', 0) >= exon_aa_start:
                label = d.get('display_name') or d.get('description', d.get('id', 'Unknown Domain'))
                overlapping_domain_names.append(label)
                overlapping_domain_details.append({
                    "label": label,
                    "name": d.get('description') or d.get('id') or d.get('hit_id') or "Domain",
                    "source": d.get('source') or d.get('type'),
                    "id": d.get('id') or d.get('hit_id'),
                    "url": d.get('url'),
                    "start": d.get('start'),
                    "end": d.get('end')
                })

    domain_count = len(overlapping_domain_names)
    cond5_no_domain = domain_count == 0
    cond6_missense = counts['missense'] < 3 + counts['frameshift'] + counts['nonsense']
    cond7_splice = counts['splice'] == 0
    cond8_no_inframe_del = counts['inframe_del'] == 0
    cond9_benign_splice = counts['benign_splice'] > 0

    # --- Step 2.5: Apply Overrides to get FINAL states ---
    # Start with original values, then overwrite if an override exists
    cond1_inframe = overrides.get('Is In-Frame', original_conds['cond1_inframe'])
    cond2_no_stop = overrides.get('No New Stop Codon', original_conds['cond2_no_stop'])
    cond3_not_terminal = overrides.get('Not First/Last Exon', original_conds['cond3_not_terminal'])
    cond4_small = overrides.get('Is <10% of Protein', original_conds['cond4_small'])
    cond5_no_domain = overrides.get('No Domain Overlap', original_conds['cond5_no_domain'])
    cond6_missense = overrides.get('Low Missense Count', original_conds['cond6_missense'])
    cond7_splice = overrides.get('No Pathogenic Splice Variants', original_conds['cond7_splice'])
    cond8_no_inframe_del = overrides.get('No Pathogenic In-Frame Deletions', original_conds['cond8_no_inframe_del'])
    cond9_benign_splice = overrides.get('Benign splice variant found', original_conds['cond9_benign_splice'])
    cond10_no_new_codon = overrides.get('No New Codon Formation', original_conds['cond10_no_new_codon'])
    has_override = bool(overrides)
    
    # Check if any override actually changed the outcome
    is_changed = False
    if has_override:
        check_map = {
            'Is In-Frame': cond1_inframe != original_conds['cond1_inframe'], 'No New Stop Codon': cond2_no_stop != original_conds['cond2_no_stop'],
            'Not First/Last Exon': cond3_not_terminal != original_conds['cond3_not_terminal'], 'Is <10% of Protein': cond4_small != original_conds['cond4_small'],
            'No Domain Overlap': cond5_no_domain != original_conds['cond5_no_domain'], 'Low Missense Count': cond6_missense != original_conds['cond6_missense'],
            'No Pathogenic Splice Variants': cond7_splice != original_conds['cond7_splice'], 'No Pathogenic In-Frame Deletions': cond8_no_inframe_del != original_conds['cond8_no_inframe_del'],
            'Benign splice variant found': cond9_benign_splice != original_conds['cond9_benign_splice'],
            'No New Codon Formation': cond10_no_new_codon != original_conds['cond10_no_new_codon']
        }
        is_changed = any(check_map.get(check_name, False) for check_name in overrides)

    # --- Step 3: Classification Logic Chain ---
    if not cond3_not_terminal:
        classification, reason = "Not Eligible", "Exon is the first or last coding exon."
    else:
        # Check consecutive exons if originally out of frame, regardless of override status of 'Is In-Frame'
        if not original_conds['cond1_inframe']:
            consecutive_res = assess_consecutive_exons(client, transcript, all_exons, target_exon, vep_entry, refseq_id_for_viewer, overrides)
            if consecutive_res: return consecutive_res

        if not cond1_inframe:
            classification, reason = "Not Eligible", "Exon is out-of-frame, which would disrupt the reading frame. Multi-exon skipping was not found to be a viable in-frame alternative."
        elif not cond2_no_stop:
            classification, reason = "Not Eligible", "Skipping this exon is predicted to create a premature stop codon."
        elif not cond7_splice:
            classification, reason = "Not Eligible", f"Exon contains {counts['splice']} pathogenic splice variant(s), indicating exon loss is pathogenic."
        elif not cond4_small and domain_count > 1:
            classification, reason = "Not Eligible", f"Exon overlaps with {domain_count} protein domains and constitutes >10% of the coding region."
        elif not cond8_no_inframe_del:
            classification, reason = "Not Eligible", f"Exon contains {counts['inframe_del']} pathogenic in-frame deletion(s)."
        elif cond9_benign_splice:
            if counts['splice'] > 0 or counts['inframe_del'] > 0:
                classification, reason = "Unlikely Eligible", "Exon has benign splice variants but also pathogenic splice or in-frame deletion variants, suggesting caution."
            else:
                classification, reason = "Eligible", "Exon contains benign splice variants, suggesting it may be safely skipped."
        elif not cond10_no_new_codon:
            classification, reason = "Unlikely Eligible", "Skipping this exon creates a novel codon at the junction that codes for a different amino acid."
        elif not cond6_missense:
            classification, reason = "Unlikely Eligible", f"Exon is a mutational hotspot with {counts['missense']} pathogenic missense variants."
        elif not cond4_small:
            classification, reason = "Unlikely Eligible", "Exon constitutes >=10% of the protein, risking major functional loss."
        elif not cond5_no_domain:
            classification, reason = "Unlikely Eligible", f"Exon overlaps with {domain_count} protein domain(s)."
        else:
            classification, reason = "Likely Eligible", "Exon meets the primary criteria for a skippable exon."
    
    note_parts = []
    if is_changed:
        note_parts.append("Classification based on manual override of guideline checks.")
    exon_len_nt = target_exon['end'] - target_exon['start'] + 1
    if exon_len_nt > 300:
        note_parts.append("Exon length exceeds 300 nucleotides. Proceed with caution.")
    note_parts.append("Please note that some exceptions for exon skipping apply. Please refer to the methods section for more details.")
    
    note = " ".join(note_parts) if note_parts else None

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
                            "name": domain.get('display_name', domain.get('description', domain.get('id', 'Domain')))
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
        "exon_length": exon_len_nt,
        "note": note,
        "frac_cds": (
            f"{((Decimal(str(exon_aa_len)) / Decimal(str(total_protein_len))) * 100):.2f}%"
            if total_protein_len > 0
            else "N/A"
        ),
        "pathogenic_variant_counts": counts, 
        "domain_count": domain_count,
        "domain_names": overlapping_domain_names,
        "domain_details": overlapping_domain_details,
        "splice_variant_links": splice_variant_links,
        "coding_exon_number": coding_exon_number,
        "total_exon_number": target_exon['total_exon_number'], 
        "gene_id": gene_id,
        "transcript_id": transcript_id,
        "clinvar_url": clinvar_url,
        "checks": {
            "Benign splice variant found": {"passed": cond9_benign_splice, "original_passed": original_conds['cond9_benign_splice'], "overridden": 'Benign splice variant found' in overrides},
            "Is In-Frame": {"passed": cond1_inframe, "original_passed": original_conds['cond1_inframe'], "overridden": 'Is In-Frame' in overrides},
            "No New Stop Codon": {"passed": cond2_no_stop, "original_passed": original_conds['cond2_no_stop'], "overridden": 'No New Stop Codon' in overrides},
            "Not First/Last Exon": {"passed": cond3_not_terminal, "original_passed": original_conds['cond3_not_terminal'], "overridden": 'Not First/Last Exon' in overrides},
            "No Pathogenic Splice Variants": {"passed": cond7_splice, "original_passed": original_conds['cond7_splice'], "overridden": 'No Pathogenic Splice Variants' in overrides},
            "No Pathogenic In-Frame Deletions": {"passed": cond8_no_inframe_del, "original_passed": original_conds['cond8_no_inframe_del'], "overridden": 'No Pathogenic In-Frame Deletions' in overrides},
            "No Domain Overlap": {"passed": cond5_no_domain, "original_passed": original_conds['cond5_no_domain'], "overridden": 'No Domain Overlap' in overrides},
            "No New Codon Formation": {"passed": cond10_no_new_codon, "original_passed": original_conds['cond10_no_new_codon'], "overridden": 'No New Codon Formation' in overrides},
            "Low Missense Count": {"passed": cond6_missense, "original_passed": original_conds['cond6_missense'], "overridden": 'Low Missense Count' in overrides},
            "Is <10% of Protein": {"passed": cond4_small, "original_passed": original_conds['cond4_small'], "overridden": 'Is <10% of Protein' in overrides}
        },
        "visualization": visualization_data
    }

def process_single_variant(query: str, client: EnsemblClient, splice_user_input: Optional[str] = None, moa_user_input: Optional[str] = None, exonic_pathogenic_user_input: Optional[str] = None, exon_skipping_overrides: Dict = {}, knockdown_overrides: Dict = {}, wt_upregulation_overrides: Dict = {}, gene_characteristics_overrides: Optional[Dict] = None) -> Dict[str, Any]:
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
        if gene_characteristics_overrides:
            gene_characteristics.update(gene_characteristics_overrides)

        protein_effect = None

        # Attempt 1: Check the chosen target_consequence directly
        if 'hgvsp' in target_consequence:
            protein_effect = target_consequence['hgvsp'].split(':')[-1]
        
        # Attempt 2: If missing, look at ALL consequences for this gene to find *any* valid HGVSp
        # (This handles cases where 'choose_best_consequence' picked a transcript that lacks p. notation)
        if not protein_effect:
            for cons in all_consequences:
                if cons.get('gene_symbol') == gene_symbol and 'hgvsp' in cons:
                    protein_effect = cons['hgvsp'].split(':')[-1]
                    break
        
        # Attempt 3: Construct from amino_acids/codons if HGVSp is still missing
        if not protein_effect:
            aa = target_consequence.get('amino_acids')
            pos = target_consequence.get('protein_start')
            if aa and pos:
                protein_effect = f"p.({pos}{aa})"

        # Attempt 4: Fallback to Consequence Type
        if not protein_effect:
            terms = target_consequence.get('consequence_terms', [])
            if terms:
                protein_effect = ", ".join(terms).replace('_', ' ').title()
            else:
                protein_effect = "Unable to calculate"

        if protein_effect:
            protein_effect = protein_effect.replace('(', '').replace(')', '').replace('/', '>')

        final_result = {
            "summary": {
                "gene": gene_symbol, 
                "transcript_id": definitive_transcript_id,
                "protein_effect": protein_effect, 
                **gene_characteristics
            },
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
        resolved_moa = moa_user_input if moa_user_input in ("GoF", "LoF", "DN") else None
        # Expose resolved_moa for frontend rendering while preserving original Known Mechanism
        final_result["summary"]["resolved_moa"] = resolved_moa
        if resolved_moa:
            moi = set(gene_characteristics.get("moi", []))
            
            if resolved_moa in("GoF", "DN"):
                is_ad = any("autosomal dominant" or "x-linked inheritance" in str(m).lower() for m in moi)
                if is_ad:
                    n1c_gene_kd = n1c_gene_knockdown_entry(gene_symbol)
                    if n1c_gene_kd:
                        final_result["assessments"]["N1C_Gene_Knockdown"] = n1c_gene_kd
                        return final_result  
                                        
                    # Only run regular knockdown if no N1C gene-level knockdown entry
                    knockdown_assessment = assess_knockdown(gene_characteristics, overrides=knockdown_overrides)
                    if resolved_moa == "DN":
                        knockdown_assessment["note"] = "Dominant Negative mechanism; allele-specific knockdown is necessary."

                    if resolved_moa == "GoF":
                        knockdown_assessment["note"] = "autosomal dominant mode of inheritance; allele-specific knockdown is recommended."
                    final_result["assessments"]["Knockdown"] = knockdown_assessment

            if resolved_moa in("GoF"):

                is_xlinked = any("x-linked inheritance" or "autosomal recessive" in str(m).lower() for m in moi)
                if is_xlinked:
                    n1c_gene_kd = n1c_gene_knockdown_entry(gene_symbol)
                    if n1c_gene_kd:
                        final_result["assessments"]["N1C_Gene_Knockdown"] = n1c_gene_kd
                        return final_result
                    
                    knockdown_assessment = assess_knockdown(gene_characteristics, overrides=knockdown_overrides)
                    final_result["assessments"]["Knockdown"] = knockdown_assessment
            

            
            elif resolved_moa == "LoF":
                # Show WT Upregulation only when LoF and Autosomal Dominant MOI (case-insensitive match)
                is_lof_ad = any("autosomal dominant inheritance" in str(m).lower() for m in moi)
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
                details = {"Confirmation Method": "User-provided validation (qPCR/RNA-seq/cDNA)"}
                splice_assessment = _evaluate_splice_variant_position(variant_identifier_from_vep, vep_entry, details, exonic_pathogenic_user_input=exonic_pathogenic_user_input)
            elif splice_user_input == 'no':
                splice_assessment = {"classification": "Unable to Assess", "reason": "User confirmed no known splice-altering effect."}
            else:
                splice_assessment = assess_splice_switching(variant_identifier_from_vep, vep_entry, gene_symbol, exonic_pathogenic_user_input=exonic_pathogenic_user_input)
            
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
                    exon_skip_result = assess_single_exon(client, query, transcript_data, all_exons, target_exon, vep_entry, overrides=exon_skipping_overrides, refseq_id_for_viewer=refseq_id_for_viewer)
                    if "visualization" in exon_skip_result and exon_skip_result["visualization"]:
                        final_result["visualization"] = exon_skip_result.pop("visualization")
                    # N1C registry exon-skipping support: if N1C lists exon skipping for this exon, mark eligible
                    try:
                        n1c_exons, n1c_links, n1c_exon_link_map = n1c_exon_skipping_exon_numbers_for_gene(gene_symbol)
                    except Exception:
                        n1c_exons, n1c_links, n1c_exon_link_map = set(), [], {}
                    try:
                        n1c_variant_exon_map = n1c_exon_skipping_variant_exon_map(gene_symbol, all_exons, client)
                    except Exception:
                        n1c_variant_exon_map = {}

                    exon_num = target_exon.get('total_exon_number')
                    matched_links: List[str] = []
                    variant_links = n1c_variant_exon_map.get(exon_num, []) if exon_num is not None else []
                    has_text_match = exon_num in n1c_exons if exon_num is not None else False
                    if variant_links:
                        matched_links.extend(variant_links)
                    elif has_text_match:
                        matched_links.extend(n1c_exon_link_map.get(exon_num, []))
                        if not matched_links and n1c_links:
                            matched_links.extend(n1c_links)

                    if variant_links or has_text_match:
                        exon_skip_result = dict(exon_skip_result)
                        note = f"N1C registry contains a project that maps to exon {exon_num} in {gene_symbol}. The input variant lies in the same exon and is therefore considered skippable."
                        prev_reason = exon_skip_result.get('reason', '')
                        exon_skip_result['reason'] = (prev_reason + ' ' + note).strip()
                        exon_skip_result['classification'] = 'Eligible'
                        det = exon_skip_result.get('details', {}) if isinstance(exon_skip_result.get('details'), dict) else {}
                        if matched_links:
                            det['N1C Exon Skipping (same exon)'] = ', '.join(sorted(set(matched_links)))
                        elif n1c_links:
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
    exonic_pathogenic_user_input = data.get('exonic_pathogenic_user_input', None)
    moa_input = data.get('moa_user_input', None)
    es_overrides = data.get('exon_skipping_overrides', {})
    kd_overrides = data.get('knockdown_overrides', {})
    wt_overrides = data.get('wt_upregulation_overrides', {})
    client = EnsemblClient()
    result = process_single_variant(query, client, splice_user_input=splice_input, moa_user_input=moa_input, exonic_pathogenic_user_input=exonic_pathogenic_user_input, exon_skipping_overrides=es_overrides, knockdown_overrides=kd_overrides, wt_upregulation_overrides=wt_overrides)
    
    return jsonify(result)
@app.route('/batch_assess', methods=['POST'])
def batch_assess():
    if 'file' not in request.files:
        return jsonify({"error": "No file part"}), 400
    file = request.files['file']
    if file.filename == '':
        return jsonify({"error": "No selected file"}), 400

    # Helper to normalize Y/N
    def norm_yn(val):
        if pd.isna(val): return None
        s = str(val).strip().upper()
        if s.startswith('Y'): return 'yes'
        if s.startswith('N'): return 'no'
        return None

    # Helper to normalize MoA
    def norm_moa(val):
        if pd.isna(val): return None
        s = str(val).strip().upper()
        if 'GOF' in s: return 'GoF'
        if 'LOF' in s: return 'LoF'
        if 'DN' in s: return 'DN'
        return None

    # Helper to normalize MOI
    def norm_moi(val):
        if pd.isna(val): return None
        s = str(val).strip().upper()
        res = []
        if 'AD' in s: res.append('Autosomal dominant inheritance')
        if 'AR' in s: res.append('Autosomal recessive inheritance')
        if 'XD' in s: res.append('X-linked dominant inheritance')
        if 'XR' in s: res.append('X-linked recessive inheritance')
        return res if res else None

    # Helper to normalize Haplo
    def norm_haplo(val):
        if pd.isna(val): return None
        s = str(val).strip().upper()
        if s.startswith('Y'): return {"text": "Sufficient evidence", "url": None}
        if s.startswith('N'): return {"text": "No evidence", "url": None}
        return None

    try:
        if file.filename.endswith('.xlsx'):
            df = pd.read_excel(file, header=None)
        else:
            df = pd.read_csv(file, header=None)
        
        # Check first row values to detect template
        first_row = df.iloc[0].astype(str).str.strip().tolist()
        is_template = 'HGVS_Variant' in first_row
        
        if is_template:
            df.columns = first_row
            df = df.iloc[1:].reset_index(drop=True)
            rows_to_process = df.to_dict('records')
        else:
            variants = df[0].dropna().astype(str).tolist()
            rows_to_process = [{'HGVS_Variant': v} for v in variants]

    except Exception as e:
        return jsonify({"error": f"Error reading file: {e}"}), 400

    client = EnsemblClient()
    output_rows = []

    for row_data in rows_to_process:
        variant = str(row_data.get('HGVS_Variant', '')).strip()
        if not variant or variant.lower() == 'nan': continue

        # Capture raw inputs for reporting
        raw_moa = row_data.get('Pathomechanism(GOF/LOF/DN)')
        raw_splice = row_data.get('Splicing effects (Y/N)')
        raw_exonic = row_data.get('Splicing-independent effects(Y/N)')
        raw_moi = row_data.get('Mode of Inheritance (AD/AR/XD/XR)')
        raw_haplo = row_data.get('Haploinsufficiency(Y/N)')
        
        user_inputs = []
        if not pd.isna(raw_moa): user_inputs.append(f"MoA: {raw_moa}")
        if not pd.isna(raw_splice): user_inputs.append(f"Splice: {raw_splice}")
        if not pd.isna(raw_exonic): user_inputs.append(f"Splice-Indep: {raw_exonic}")
        if not pd.isna(raw_moi): user_inputs.append(f"MOI: {raw_moi}")
        if not pd.isna(raw_haplo): user_inputs.append(f"Haplo: {raw_haplo}")
        user_input_str = "; ".join(user_inputs) if user_inputs else "None"

        # Extract overrides
        moa_input = norm_moa(raw_moa)
        splice_input = norm_yn(raw_splice)
        exonic_pathogenic_input = norm_yn(raw_exonic)
        
        moi_override = norm_moi(raw_moi)
        haplo_override = norm_haplo(raw_haplo)

        gene_char_overrides = {}
        if moi_override: gene_char_overrides['moi'] = moi_override
        if haplo_override: gene_char_overrides['haploinsufficiency'] = haplo_override

        result = process_single_variant(
            variant, client, 
            splice_user_input=splice_input, 
            moa_user_input=moa_input, 
            exonic_pathogenic_user_input=exonic_pathogenic_input,
            gene_characteristics_overrides=gene_char_overrides
        )
        
        row = {"Variant": variant}
        row["User Input"] = user_input_str
        summary = result.get("summary", {})
        assessments = result.get("assessments", {})

        row["Gene"] = summary.get("gene", "N/A")
        row["ClinGen MOI"] = ', '.join(summary.get("moi", []))
        row["ClinGen MOA"] = ', '.join(summary.get("moa", []))
        # Ensembl Transcript ID and link
        transcript_id = summary.get("transcript_id")
        row["Ensembl Transcript"] = transcript_id or "N/A"
        row["Ensembl Transcript Link"] = (
            f"https://www.ensembl.org/Homo_sapiens/Transcript/Summary?t={transcript_id}" if transcript_id else "N/A"
        )
        
        haplo_info = summary.get("haploinsufficiency", {})
        row["Haploinsufficiency"] = haplo_info.get("text", "N/A")
        row["ClinGen Link"] = haplo_info.get("url", "N/A")
        
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
            if isinstance(status, dict):
                row[f"ES Check: {check}"] = status.get("original_passed", "N/A")
            else:
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
            domain_details = skip.get("domain_details") or []
            if domain_details:
                formatted_domains = []
                for d in domain_details:
                    label = d.get("label") or d.get("name") or "Domain"
                    url = d.get("url")
                    formatted_domains.append(f"{label} [{url}]" if url else label)
                row["Domains"] = "; ".join(formatted_domains)
            else:
                row["Domains"] = ", ".join(domain_names) if domain_names else "N/A"

        # Splice Correction
        splice = assessments.get("Splice_Switching", {})
        splice_class = splice.get("classification", "Unable to Assess")
        if splice_class == "Awaiting User Input":
            splice_class = "Splicing-independent effects unknown. Manual/Single Variant processing needed"
        row["Splice Correction Assessment"] = splice_class
        row["Splicing Validation DOI"] = splice.get("details", {}).get("Publication DOI", "NA")
        # Splicing DB/Source links if available
        splice_details = splice.get("details", {}) if isinstance(splice.get("details", {}), dict) else {}
        splicing_db_link = splice_details.get("SSCVDB Gene Page") or splice_details.get("Publication") or "N/A"
        row["Splicing DB Link"] = splicing_db_link

        # WT Upregulation and Knockdown (assessments remain)
        row["WT-Upregulation"] = wt_up.get("classification", "NA")
        row["Knockdown"] = assessments.get("Knockdown", {}).get("classification", "NA")

        # Manual validation needs and dual MoA assessment when MoA unresolved
        manual_needs = []
        summary_moa_list = summary.get("moa", []) or []
        resolved_moa = summary.get("resolved_moa")
        # Splice manual validation prompt
        if splice.get("user_validation_prompt") or (splice.get("classification") in ("Not in Database", "Unable to Assess")):
            manual_needs.append("Splice validation (Variant was not found in SpliceVarDB/SSCVDB and therefore requires user confirmation)")
        # MoA unclear -> warn
        if not resolved_moa and (len(summary_moa_list) != 1):
            manual_needs.append("Mechanism selection (DN vs GoF vs LoF)")

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
    
@app.route('/download_batch_template')
def download_batch_template():
    """Serves the batch processing template file."""
    try:
        BASE_DIR = os.path.dirname(os.path.abspath(__file__))
        DATA_DIR = os.path.join(BASE_DIR, 'data')
        template_path = os.path.join(DATA_DIR, 'batch_template.xlsx')
        return send_file(
            template_path,
            as_attachment=True,
            download_name='AVEC_batch_template.xlsx',
            mimetype='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'
        )
    except FileNotFoundError:
        return "Template file not found on server.", 404



