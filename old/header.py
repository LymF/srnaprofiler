#!/usr/bin/env python3
"""
Alternative header designs for sRNA Profiler
Choose your favorite and replace the print_header() function
"""

def print_header_v1():
    """Version 1: Compact with colors matching logo style"""
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    WHITE = '\033[97m'
    BOLD = '\033[1m'
    RESET = '\033[0m'
    BG_BLUE = '\033[44m'
    BG_GREEN = '\033[42m'
    
    header = f"""
{BLUE}{'─' * 70}{RESET}
{BOLD}{BG_BLUE}{WHITE}                                                                      {RESET}
{BOLD}{BG_BLUE}{WHITE}     sRNA{RESET}{BOLD}{BG_GREEN}{WHITE}Profiler{RESET}{BOLD}{BG_BLUE}{WHITE}     🧬  Small RNA Analysis Pipeline        {RESET}
{BOLD}{BG_BLUE}{WHITE}                                                                      {RESET}
{BLUE}{'─' * 70}{RESET}

{GREEN}  ╭─────────────────────────────────────────────────────────────╮
  │ {WHITE}📊 Quality Control   {GREEN}│{WHITE} 🧭 Read Mapping   {GREEN}│{WHITE} 🔬 Cluster Analysis{GREEN} │
  ╰─────────────────────────────────────────────────────────────╯{RESET}

{BLUE}  Version 1.0{RESET}  |  {GREEN}Bioinformatics Tool for sRNA-Seq Analysis{RESET}
{BLUE}{'─' * 70}{RESET}
"""
    print(header)


def print_header_v2():
    """Version 2: Large ASCII art logo"""
    BLUE = '\033[34m'
    GREEN = '\033[92m'
    DARK_GREEN = '\033[32m'
    WHITE = '\033[97m'
    BOLD = '\033[1m'
    RESET = '\033[0m'
    
    header = f"""
{BOLD}{BLUE}╔═══════════════════════════════════════════════════════════════════════╗{RESET}
{BOLD}{BLUE}║                                                                       ║{RESET}
{BOLD}{BLUE}║  {WHITE}     ____  _   _    _    {GREEN}____            __ _ _           {BLUE}         ║{RESET}
{BOLD}{BLUE}║  {WHITE} ___||  _ \\| \\ | |  / \\  {GREEN}|  _ \\ _ __ ___ / _(_) | ___ _ __ {BLUE}         ║{RESET}
{BOLD}{BLUE}║  {WHITE}/ __|| |_) |  \\| | / _ \\ {GREEN}| |_) | '__/ _ \\| |_| | |/ _ \\ '__|{BLUE}         ║{RESET}
{BOLD}{BLUE}║  {WHITE}\\__ \\|  _ <| |\\  |/ ___ \\{GREEN}|  __/| | | (_) |  _| | |  __/ |  {BLUE}         ║{RESET}
{BOLD}{BLUE}║  {WHITE}|___/|_| \\_\\_| \\_/_/   \\_\\{GREEN}_|   |_|  \\___/|_| |_|_|\\___|_|  {BLUE}         ║{RESET}
{BOLD}{BLUE}║                                                                       ║{RESET}
{BOLD}{BLUE}║  {WHITE}🧬  Small RNA Sequencing Analysis Pipeline  🧬{BLUE}                      ║{RESET}
{BOLD}{BLUE}╚═══════════════════════════════════════════════════════════════════════╝{RESET}
{GREEN}         QC & Trimming  •  Read Mapping  •  Cluster Analysis{RESET}
"""
    print(header)


def print_header_v3():
    """Version 3: Minimalist with pill-shaped design"""
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    WHITE = '\033[97m'
    CYAN = '\033[96m'
    BOLD = '\033[1m'
    RESET = '\033[0m'
    
    header = f"""
    {BLUE}╭{'─' * 50}╮{RESET}
    {BLUE}│{BOLD}{WHITE}  sRNA{GREEN}Profiler{WHITE}  🧬  Small RNA Analysis Pipeline {RESET}{BLUE}│{RESET}
    {BLUE}╰{'─' * 50}╯{RESET}
    
    {CYAN}Quality Control │ Read Mapping │ Cluster Analysis{RESET}
    {BLUE}{'━' * 52}{RESET}
"""
    print(header)


def print_header_v4():
    """Version 4: Box drawing with gradient effect"""
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    DARK_GREEN = '\033[32m'
    WHITE = '\033[97m'
    BOLD = '\033[1m'
    DIM = '\033[2m'
    RESET = '\033[0m'
    
    header = f"""
{BLUE}┏━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃                                                                    ┃
┃   {BOLD}{WHITE}███████╗{BLUE}██████╗ {WHITE}███╗   ██╗{GREEN} █████╗ {DARK_GREEN}██████╗ ██████╗  ██████╗ {BLUE}   ┃
┃   {BOLD}{WHITE}██╔════╝{BLUE}██╔══██╗{WHITE}████╗  ██║{GREEN}██╔══██╗{DARK_GREEN}██╔══██╗██╔══██╗██╔═══██╗{BLUE}   ┃
┃   {BOLD}{WHITE}███████╗{BLUE}██████╔╝{WHITE}██╔██╗ ██║{GREEN}███████║{DARK_GREEN}██████╔╝██████╔╝██║   ██║{BLUE}   ┃
┃   {BOLD}{WHITE}╚════██║{BLUE}██╔══██╗{WHITE}██║╚██╗██║{GREEN}██╔══██║{DARK_GREEN}██╔═══╝ ██╔══██╗██║   ██║{BLUE}   ┃
┃   {BOLD}{WHITE}███████║{BLUE}██║  ██║{WHITE}██║ ╚████║{GREEN}██║  ██║{DARK_GREEN}██║     ██║  ██║╚██████╔╝{BLUE}   ┃
┃   {BOLD}{WHITE}╚══════╝{BLUE}╚═╝  ╚═╝{WHITE}╚═╝  ╚═══╝{GREEN}╚═╝  ╚═╝{DARK_GREEN}╚═╝     ╚═╝  ╚═╝ ╚═════╝ {BLUE}   ┃
┃                                                                    ┃
┃              {DIM}{WHITE}Small RNA Sequencing Analysis Pipeline{RESET}{BLUE}               ┃
┗━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┛{RESET}
"""
    print(header)


def print_header_v5():
    """Version 5: Clean and professional with DNA helix symbol"""
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    WHITE = '\033[97m'
    BOLD = '\033[1m'
    RESET = '\033[0m'
    BG_BLUE = '\033[104m'  # Light blue background
    
    header = f"""
{BOLD}{BLUE}                    ╔══════════════════════════════╗
                    ║  {WHITE}sRNA{GREEN}Profiler{BLUE}  │  v1.0       ║
                    ╚══════════════════════════════╝{RESET}

                {WHITE}🧬  Small RNA Sequencing Analysis Pipeline  🧬{RESET}

{BLUE}    ┌─────────────────────────────────────────────────────────┐
    │ {WHITE}• Quality Control & Read Trimming                        {BLUE}│
    │ {WHITE}• Genome Preparation & Tiling                            {BLUE}│
    │ {WHITE}• Read Alignment with Bowtie                             {BLUE}│
    │ {WHITE}• Cluster Identification & Classification                {BLUE}│
    │ {WHITE}• Pattern Analysis & Visualization                       {BLUE}│
    └─────────────────────────────────────────────────────────┘{RESET}
"""
    print(header)


def print_header_v6():
    """Version 6: Rounded corners matching the logo style"""
    BLUE = '\033[34m'
    GREEN = '\033[92m'
    WHITE = '\033[97m'
    BOLD = '\033[1m'
    RESET = '\033[0m'
    BG_BLUE = '\033[44m'
    BG_GREEN = '\033[42m'
    
    # Using Unicode block characters for a pill-shaped effect
    header = f"""
    {BOLD}{BLUE}╭{'═' * 46}╮{RESET}
    {BOLD}{BLUE}║{BG_BLUE}{WHITE}                                              {RESET}{BOLD}{BLUE}║{RESET}
    {BOLD}{BLUE}║{BG_BLUE}{WHITE}     sRNA{RESET}{BOLD}{BG_GREEN}{WHITE}Profiler{RESET}{BOLD}{BG_BLUE}{WHITE}   🧬                     {RESET}{BOLD}{BLUE}║{RESET}
    {BOLD}{BLUE}║{BG_BLUE}{WHITE}     Small RNA Analysis Pipeline              {RESET}{BOLD}{BLUE}║{RESET}
    {BOLD}{BLUE}║{BG_BLUE}{WHITE}                                              {RESET}{BOLD}{BLUE}║{RESET}
    {BOLD}{BLUE}╰{'═' * 46}╯{RESET}
    
    {GREEN}  QC & Trimming  •  Read Mapping  •  Cluster Analysis{RESET}
"""
    print(header)


def print_header_v7():
    """Version 7: Simple with emoji and clean lines"""
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    WHITE = '\033[97m'
    BOLD = '\033[1m'
    ITALIC = '\033[3m'
    RESET = '\033[0m'
    
    header = f"""
{BLUE}{'━' * 60}{RESET}

       {BOLD}{BLUE}sRNA{GREEN}Profiler{RESET}  🧬  {ITALIC}Small RNA Analysis Pipeline{RESET}

{BLUE}{'━' * 60}{RESET}

   📊 QC & Trimming    🧭 Mapping    🔬 Cluster Analysis

{BLUE}{'━' * 60}{RESET}
"""
    print(header)


# Test all versions
if __name__ == "__main__":
    print("\n\n=== VERSION 1 ===")
    print_header_v1()
    
    print("\n\n=== VERSION 2 ===")
    print_header_v2()
    
    print("\n\n=== VERSION 3 ===")
    print_header_v3()
    
    print("\n\n=== VERSION 4 ===")
    print_header_v4()
    
    print("\n\n=== VERSION 5 ===")
    print_header_v5()
    
    print("\n\n=== VERSION 6 ===")
    print_header_v6()
    
    print("\n\n=== VERSION 7 ===")
    print_header_v7()
    
    print("\n\nEscolha sua versão favorita e substitua a função print_header() no srnaprofiler.py!")