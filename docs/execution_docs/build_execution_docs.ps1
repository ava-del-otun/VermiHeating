Set-StrictMode -Version Latest
$ErrorActionPreference = 'Stop'

function Resolve-PdfLatex {
    $command = Get-Command 'pdflatex' -ErrorAction SilentlyContinue
    if ($command) {
        return $command.Source
    }

    $candidates = @(
        'C:\Users\cezar\AppData\Local\Programs\MiKTeX\miktex\bin\x64\pdflatex.exe',
        'C:\Program Files\MiKTeX\miktex\bin\x64\pdflatex.exe'
    )

    foreach ($candidate in $candidates) {
        if (Test-Path $candidate) {
            return $candidate
        }
    }

    throw 'Could not locate pdflatex.exe. Install MiKTeX or add pdflatex to PATH.'
}

$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$pdflatex = Resolve-PdfLatex
$documents = @(
    (Join-Path $scriptDir 'physics\physics_and_math_reference.tex'),
    (Join-Path $scriptDir 'physics\implementation_traceability.tex')
)

foreach ($document in $documents) {
    $documentDir = Split-Path -Parent $document
    $documentName = Split-Path -Leaf $document
    $baseName = [System.IO.Path]::GetFileNameWithoutExtension($documentName)
    $logPath = Join-Path $documentDir ($baseName + '.log')
    $pdfPath = Join-Path $documentDir ($baseName + '.pdf')

    Push-Location $documentDir
    try {
        foreach ($pass in 1..2) {
            $logTimestampBefore = if (Test-Path $logPath) { (Get-Item $logPath).LastWriteTimeUtc } else { $null }
            $pdfTimestampBefore = if (Test-Path $pdfPath) { (Get-Item $pdfPath).LastWriteTimeUtc } else { $null }

            Write-Host "Compiling $documentName (pass $pass of 2)..."
            & $pdflatex '-interaction=nonstopmode' '-halt-on-error' $documentName
            if (-not (Test-Path $logPath)) {
                throw "pdflatex did not produce $logPath."
            }

            $logTimestampAfter = (Get-Item $logPath).LastWriteTimeUtc
            $pdfTimestampAfter = if (Test-Path $pdfPath) { (Get-Item $pdfPath).LastWriteTimeUtc } else { $null }
            $logHasOutput = Select-String -Path $logPath -Pattern 'Output written on ' -Quiet
            $logHasFatal = Select-String -Path $logPath -Pattern 'Fatal error occurred|Emergency stop' -Quiet
            $logUpdated = ($null -eq $logTimestampBefore) -or ($logTimestampAfter -gt $logTimestampBefore)
            $pdfUpdated = ($null -ne $pdfTimestampAfter) -and (($null -eq $pdfTimestampBefore) -or ($pdfTimestampAfter -gt $pdfTimestampBefore))

            if (-not $logHasOutput -or $logHasFatal -or -not ($logUpdated -or $pdfUpdated)) {
                throw "pdflatex failed for $documentName on pass $pass."
            }
        }
    }
    finally {
        Pop-Location
    }
}
